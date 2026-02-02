"""Resilience capacity of burning areas using Landsat Collection 2 time series.

Analyze vegetation recovery in fire-affected areas by computing annual
Normalized Burn Ratio (NBR) and Normalized Difference Vegetation Index (NDVI)
time series from Landsat Collection 2 Level-2 surface reflectance data.

The workflow merges imagery from Landsat 5 TM, 7 ETM+, 8 OLI, and 9 OLI-2,
applies Collection 2 scale factors and QA_PIXEL cloud/shadow masking, then
produces per-year median composites within a configurable seasonal window
(default July--August, roughly corresponding to peak vegetation expansion in
the Mediterranean and temperate zones).

Recovery is quantified per fire site by comparing pre-fire baseline index
values with post-fire trajectories, yielding a dimensionless recovery ratio.

Examples
--------
Run the built-in demo (Dollar Lake fire, Mt. Hood, OR, 2011)::

    $ uv run res2fire.py --project MY_GEE_PROJECT

Process a shapefile of fire perimeters::

    $ uv run res2fire.py -f fires.shp -o results.csv --project MY_GEE_PROJECT

Notes
-----
* Requires an authenticated Google Earth Engine account.
  Run ``earthengine authenticate`` before first use.
* Landsat Collection 2 surface reflectance data are already cross-calibrated
  across sensors, so no manual ETM+-to-OLI harmonization is needed.

References
----------
.. [1] USGS, "Landsat Collection 2 Level-2 Science Products",
       https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products
.. [2] Key, C.H. and Benson, N.C. (2006), "Landscape Assessment (LA)",
       FIREMON: Fire Effects Monitoring and Inventory System, RMRS-GTR-164-CD.
"""

import argparse
import csv
import sys

import ee
import geopandas as gpd
from shapely.geometry import mapping

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LANDSAT_COLLECTIONS: dict[str, str] = {
    "L5": "LANDSAT/LT05/C02/T1_L2",
    "L7": "LANDSAT/LE07/C02/T1_L2",
    "L8": "LANDSAT/LC08/C02/T1_L2",
    "L9": "LANDSAT/LC09/C02/T1_L2",
}
"""Mapping of sensor short-names to GEE Collection 2 Tier-1 Level-2 asset IDs."""

OLI_BANDS: list[str] = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"]
"""OLI / OLI-2 surface reflectance band names (Landsat 8 & 9)."""

ETM_BANDS: list[str] = ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"]
"""TM / ETM+ surface reflectance band names (Landsat 5 & 7)."""

STD_NAMES: list[str] = ["Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"]
"""Standardized wavelength-based band names used after renaming."""


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

def apply_scale_factors(image: ee.Image) -> ee.Image:
    """Apply Landsat Collection 2 surface reflectance scale factors.

    Converts raw integer DN values to physical reflectance in the range
    [0, 1] using the USGS-provided linear transformation::

        reflectance = DN * 0.0000275 - 0.2

    Parameters
    ----------
    image : ee.Image
        A raw Landsat Collection 2 Level-2 image.

    Returns
    -------
    ee.Image
        The same image with optical ``SR_B*`` bands rescaled to reflectance.
    """
    optical = image.select("SR_B.").multiply(0.0000275).add(-0.2)
    return image.addBands(optical, overwrite=True)


def cloud_mask_c2(image: ee.Image) -> ee.Image:
    """Mask clouds and cloud shadows using the Collection 2 QA_PIXEL band.

    The following QA_PIXEL bit flags are used to identify and mask
    contaminated pixels:

    * Bit 1 -- Dilated cloud
    * Bit 2 -- Cirrus
    * Bit 3 -- Cloud
    * Bit 4 -- Cloud shadow

    Parameters
    ----------
    image : ee.Image
        A Landsat Collection 2 image containing a ``QA_PIXEL`` band.

    Returns
    -------
    ee.Image
        The input image with cloud/shadow pixels masked out.
    """
    qa = image.select("QA_PIXEL")
    dilated_cloud = 1 << 1
    cirrus = 1 << 2
    cloud = 1 << 3
    cloud_shadow = 1 << 4
    mask = (
        qa.bitwiseAnd(dilated_cloud).eq(0)
        .And(qa.bitwiseAnd(cirrus).eq(0))
        .And(qa.bitwiseAnd(cloud).eq(0))
        .And(qa.bitwiseAnd(cloud_shadow).eq(0))
    )
    return image.updateMask(mask)


def rename_oli(image: ee.Image) -> ee.Image:
    """Rename Landsat 8/9 OLI bands to standardized wavelength names.

    Parameters
    ----------
    image : ee.Image
        An OLI or OLI-2 image with original ``SR_B*`` band names.

    Returns
    -------
    ee.Image
        Image with bands renamed to :data:`STD_NAMES`.
    """
    return image.select(OLI_BANDS, STD_NAMES)


def rename_etm(image: ee.Image) -> ee.Image:
    """Rename Landsat 5/7 TM/ETM+ bands to standardized wavelength names.

    Parameters
    ----------
    image : ee.Image
        A TM or ETM+ image with original ``SR_B*`` band names.

    Returns
    -------
    ee.Image
        Image with bands renamed to :data:`STD_NAMES`.
    """
    return image.select(ETM_BANDS, STD_NAMES)


def add_indices(image: ee.Image) -> ee.Image:
    """Compute and append NBR and NDVI spectral index bands.

    Indices are calculated as normalized differences:

    * **NBR** = (NIR - SWIR2) / (NIR + SWIR2)  [1]_
    * **NDVI** = (NIR - Red) / (NIR + Red)

    Parameters
    ----------
    image : ee.Image
        An image with standardized band names (``NIR``, ``SWIR2``, ``Red``).

    Returns
    -------
    ee.Image
        The input image with ``NBR`` and ``NDVI`` bands appended.

    References
    ----------
    .. [1] Key, C.H. and Benson, N.C. (2006), FIREMON: Fire Effects
           Monitoring and Inventory System, RMRS-GTR-164-CD.
    """
    nbr = image.normalizedDifference(["NIR", "SWIR2"]).rename("NBR")
    ndvi = image.normalizedDifference(["NIR", "Red"]).rename("NDVI")
    return image.addBands([nbr, ndvi])


def prep_oli(image: ee.Image) -> ee.Image:
    """Full preprocessing pipeline for Landsat 8/9 OLI imagery.

    Applies, in order: cloud masking, scale factor correction, band
    renaming, and spectral index computation.

    Parameters
    ----------
    image : ee.Image
        A raw Landsat 8 or 9 Collection 2 Level-2 image.

    Returns
    -------
    ee.Image
        Cloud-masked image with reflectance values and NBR/NDVI bands.
    """
    return add_indices(rename_oli(apply_scale_factors(cloud_mask_c2(image))))


def prep_etm(image: ee.Image) -> ee.Image:
    """Full preprocessing pipeline for Landsat 5/7 TM/ETM+ imagery.

    Applies, in order: cloud masking, scale factor correction, band
    renaming, and spectral index computation.

    Parameters
    ----------
    image : ee.Image
        A raw Landsat 5 or 7 Collection 2 Level-2 image.

    Returns
    -------
    ee.Image
        Cloud-masked image with reflectance values and NBR/NDVI bands.
    """
    return add_indices(rename_etm(apply_scale_factors(cloud_mask_c2(image))))


# ---------------------------------------------------------------------------
# Time-series construction
# ---------------------------------------------------------------------------

def get_merged_collection(
    geometry: ee.Geometry,
    start_year: int,
    end_year: int,
    month_start: int = 7,
    month_end: int = 8,
) -> ee.ImageCollection:
    """Build a merged, preprocessed Landsat image collection.

    Loads Landsat 5, 7, 8, and 9 Collection 2 Tier-1 Level-2 imagery,
    applies the appropriate preprocessing pipeline to each sensor, and
    merges the results into a single collection filtered by spatial
    extent, date range, and seasonal window.

    Parameters
    ----------
    geometry : ee.Geometry
        Area of interest used for spatial filtering.
    start_year : int
        First year of the analysis period (inclusive).
    end_year : int
        Last year of the analysis period (inclusive).
    month_start : int, optional
        First month of the seasonal window (1--12). Default is ``7`` (July).
    month_end : int, optional
        Last month of the seasonal window (1--12, inclusive).
        Default is ``8`` (August).

    Returns
    -------
    ee.ImageCollection
        Merged collection with standardized band names and NBR/NDVI bands.
    """
    date_start = f"{start_year}-01-01"
    date_end = f"{end_year + 1}-01-01"

    def _load(key: str) -> ee.ImageCollection:
        col = (
            ee.ImageCollection(LANDSAT_COLLECTIONS[key])
            .filterBounds(geometry)
            .filterDate(date_start, date_end)
            .filter(ee.Filter.calendarRange(month_start, month_end, "month"))
        )
        prep = prep_oli if key in ("L8", "L9") else prep_etm
        return col.map(prep)

    merged = _load("L5").merge(_load("L7")).merge(_load("L8")).merge(_load("L9"))
    return merged


def annual_composites(
    collection: ee.ImageCollection,
    start_year: int,
    end_year: int,
    month_start: int = 7,
    month_end: int = 8,
) -> list[dict]:
    """Reduce an image collection to annual median composites.

    For each year in the range ``[start_year, end_year]``, filters the
    collection to the seasonal window and computes a pixel-wise median.
    Years with no cloud-free observations will produce an empty image
    (all pixels masked), which downstream functions handle gracefully as
    ``None`` values.

    Parameters
    ----------
    collection : ee.ImageCollection
        A preprocessed Landsat collection (output of
        :func:`get_merged_collection`).
    start_year : int
        First year (inclusive).
    end_year : int
        Last year (inclusive).
    month_start : int, optional
        Start month of the seasonal window. Default is ``7``.
    month_end : int, optional
        End month of the seasonal window (inclusive). Default is ``8``.

    Returns
    -------
    list[dict]
        A list of ``{"year": int, "image": ee.Image}`` dictionaries, one
        per year, sorted chronologically.
    """
    results = []
    for year in range(start_year, end_year + 1):
        t1 = f"{year}-{month_start:02d}-01"
        m_end = month_end + 1
        t2 = f"{year + 1}-01-01" if month_end == 12 else f"{year}-{m_end:02d}-01"
        annual = collection.filterDate(t1, t2).median()
        results.append({"year": year, "image": annual})
    return results


def extract_values_batch(
    composites: list[dict],
    geometry: ee.Geometry,
    scale: int = 30,
) -> list[dict]:
    """Extract mean NBR and NDVI from annual composites over a geometry.

    Constructs all reduction operations server-side and retrieves results
    with a single ``getInfo()`` call, avoiding the overhead of one
    round-trip per year.

    Parameters
    ----------
    composites : list[dict]
        Annual composites as returned by :func:`annual_composites`.
    geometry : ee.Geometry
        Point or polygon over which to compute the spatial mean.
    scale : int, optional
        Pixel resolution in metres for the reduction. Default is ``30``
        (native Landsat resolution).

    Returns
    -------
    list[dict]
        A list of ``{"year": int, "NBR": float | None, "NDVI": float | None}``
        dictionaries sorted by year.  ``None`` indicates no valid pixels
        were available for that year.
    """
    if not composites:
        return []

    features = []
    for comp in composites:
        image = comp["image"].select(["NBR", "NDVI"])
        stats = image.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=scale,
            maxPixels=1e9,
        )
        feat = ee.Feature(None, {"year": comp["year"]}).set(stats)
        features.append(feat)

    fc = ee.FeatureCollection(features)
    all_info = fc.getInfo()["features"]

    rows = []
    for info in all_info:
        props = info["properties"]
        rows.append({
            "year": props["year"],
            "NBR": props.get("NBR"),
            "NDVI": props.get("NDVI"),
        })
    return sorted(rows, key=lambda r: r["year"])


# ---------------------------------------------------------------------------
# Recovery metrics
# ---------------------------------------------------------------------------

def compute_recovery(
    time_series: list[dict],
    fire_year: int,
    index: str = "NBR",
) -> dict:
    """Compute vegetation recovery metrics for a single fire event.

    Metrics are derived by comparing pre-fire baseline values with the
    post-fire trajectory of a spectral vegetation index.

    The **recovery ratio** is defined as::

        RR = (V_latest - V_post_min) / (V_pre - V_post_min)

    where *V_pre* is the mean index value over the three years preceding
    the fire, *V_post_min* is the minimum post-fire value, and *V_latest*
    is the most recent valid observation.  A ratio of 1.0 indicates that
    the vegetation has returned to its pre-fire baseline.

    Parameters
    ----------
    time_series : list[dict]
        Annual index values as returned by :func:`extract_values_batch`.
    fire_year : int
        Calendar year in which the fire occurred.
    index : str, optional
        Name of the spectral index to evaluate (``"NBR"`` or ``"NDVI"``).
        Default is ``"NBR"``.

    Returns
    -------
    dict
        A dictionary with keys:

        * ``pre_fire`` (*float | None*) -- Mean index value for up to 3
          pre-fire years.
        * ``post_fire_min`` (*float | None*) -- Minimum post-fire index value.
        * ``latest`` (*float | None*) -- Index value for the most recent
          year with valid data.
        * ``latest_year`` (*int | None*) -- Year of the ``latest`` value.
        * ``delta`` (*float | None*) -- Magnitude of the fire-induced drop
          (``pre_fire - post_fire_min``).
        * ``recovery_ratio`` (*float | None*) -- Dimensionless recovery
          ratio (1.0 = full recovery).
    """
    pre = [r[index] for r in time_series if r[index] is not None and r["year"] < fire_year]
    post = {r["year"]: r[index] for r in time_series if r[index] is not None and r["year"] >= fire_year}

    if not pre or not post:
        return {"pre_fire": None, "post_fire_min": None, "latest": None,
                "latest_year": None, "delta": None, "recovery_ratio": None}

    pre_mean = sum(pre[-3:]) / len(pre[-3:])
    post_min = min(post.values())

    latest_year = max(post.keys())
    latest = post[latest_year]
    delta = pre_mean - post_min

    recovery_ratio = (latest - post_min) / delta if delta > 0 else None

    return {
        "pre_fire": round(pre_mean, 4),
        "post_fire_min": round(post_min, 4),
        "latest": round(latest, 4),
        "latest_year": latest_year,
        "delta": round(delta, 4),
        "recovery_ratio": round(recovery_ratio, 4) if recovery_ratio is not None else None,
    }


# ---------------------------------------------------------------------------
# Vegetation type summary
# ---------------------------------------------------------------------------

def _print_veg_type_summary(all_rows: list[dict]) -> None:
    """Print a summary table of recovery metrics grouped by vegetation type.

    For each vegetation type present in the results, computes the mean
    NBR and NDVI recovery ratios across all fire sites of that type and
    prints a formatted table to stdout.

    Parameters
    ----------
    all_rows : list[dict]
        The full list of output rows (one per fire site per year).
        Each row must contain ``veg_type``, ``name``, ``NBR_recovery_ratio``,
        and ``NDVI_recovery_ratio`` keys.
    """
    # Collect one recovery ratio per fire site (avoid counting each year row)
    sites: dict[str, dict[str, list[float]]] = {}
    seen: set[tuple[str, str]] = set()

    for r in all_rows:
        veg = r["veg_type"]
        name = r["name"]
        if (veg, name) in seen:
            continue
        seen.add((veg, name))

        if veg not in sites:
            sites[veg] = {"nbr": [], "ndvi": [], "count": 0}

        sites[veg]["count"] += 1
        if r.get("NBR_recovery_ratio") is not None:
            sites[veg]["nbr"].append(r["NBR_recovery_ratio"])
        if r.get("NDVI_recovery_ratio") is not None:
            sites[veg]["ndvi"].append(r["NDVI_recovery_ratio"])

    print("\n" + "=" * 64)
    print("Recovery Summary by Vegetation Type")
    print("=" * 64)
    print(f"{'Type':<20} {'Sites':>5}  {'Mean NBR RR':>11}  {'Mean NDVI RR':>12}")
    print("-" * 64)

    for veg in sorted(sites):
        s = sites[veg]
        n = s["count"]
        nbr_mean = f"{sum(s['nbr']) / len(s['nbr']):.4f}" if s["nbr"] else "N/A"
        ndvi_mean = f"{sum(s['ndvi']) / len(s['ndvi']):.4f}" if s["ndvi"] else "N/A"
        print(f"{veg:<20} {n:>5}  {nbr_mean:>11}  {ndvi_mean:>12}")

    print("=" * 64)


# ---------------------------------------------------------------------------
# Shapefile processing
# ---------------------------------------------------------------------------

def process_shapefile(
    shp_path: str,
    years_before: int,
    years_after: int,
    month_start: int,
    month_end: int,
    output_csv: str,
) -> None:
    """Process a fire-areas shapefile and write recovery results to CSV.

    Iterates over each feature in the shapefile, builds a Landsat time
    series around the fire year, extracts annual NBR/NDVI values, computes
    recovery metrics, and writes all results to a single CSV file.

    Parameters
    ----------
    shp_path : str
        Path to an OGR-readable vector file (Shapefile, GeoPackage, etc.).
        Must contain the following attribute columns:

        * ``fire_year`` (*int*, **required**) -- Year of the fire event.
        * ``veg_type`` (*str*, optional) -- Vegetation type label
          (e.g. ``"broadleaf"``, ``"conifer"``, ``"macchia"``).
        * ``name`` or ``id`` (*str*, optional) -- Human-readable feature
          identifier.

        Geometry may be Point, Polygon, or MultiPolygon.
    years_before : int
        Number of pre-fire years to include in the time series.
    years_after : int
        Number of post-fire years to include in the time series.
    month_start : int
        Start month of the seasonal compositing window (1--12).
    month_end : int
        End month of the seasonal compositing window (1--12, inclusive).
    output_csv : str
        Destination path for the output CSV file.

    Raises
    ------
    SystemExit
        If the shapefile does not contain the required ``fire_year`` column.
    """
    gdf = gpd.read_file(shp_path)

    required = {"fire_year"}
    missing = required - set(gdf.columns)
    if missing:
        print(f"Error: shapefile missing columns: {missing}")
        sys.exit(1)

    all_rows = []

    for idx, row in gdf.iterrows():
        fire_year = int(row["fire_year"])
        veg_type = row.get("veg_type", "unknown")
        name = row.get("name", row.get("id", str(idx)))
        geom = row.geometry

        start_year = fire_year - years_before
        end_year = fire_year + years_after

        print(f"Processing fire '{name}' ({veg_type}), year {fire_year}, "
              f"range {start_year}-{end_year}...")

        ee_geom = ee.Geometry(mapping(geom))
        if geom.geom_type == "Point":
            buffer = ee_geom.buffer(500)
        else:
            buffer = ee_geom

        collection = get_merged_collection(buffer, start_year, end_year,
                                           month_start, month_end)
        composites = annual_composites(collection, start_year, end_year,
                                       month_start, month_end)
        ts = extract_values_batch(composites, ee_geom)

        recovery_nbr = compute_recovery(ts, fire_year, "NBR")
        recovery_ndvi = compute_recovery(ts, fire_year, "NDVI")

        for r in ts:
            all_rows.append({
                "name": name,
                "veg_type": veg_type,
                "fire_year": fire_year,
                **r,
                **{f"NBR_{k}": v for k, v in recovery_nbr.items()},
                **{f"NDVI_{k}": v for k, v in recovery_ndvi.items()},
            })

    if all_rows:
        keys = all_rows[0].keys()
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"\nResults written to {output_csv}")
        _print_veg_type_summary(all_rows)
    else:
        print("No results produced.")


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def run_demo(
    years_before: int,
    years_after: int,
    month_start: int,
    month_end: int,
    output_csv: str = "fire_recovery.csv",
) -> None:
    """Run a demonstration analysis on the 2011 Dollar Lake fire (Mt. Hood, OR).

    Prints an annual NBR/NDVI table and recovery metrics to stdout, and
    writes results to a CSV file.

    Parameters
    ----------
    years_before : int
        Number of pre-fire years.
    years_after : int
        Number of post-fire years.
    month_start : int
        Start month of the seasonal window.
    month_end : int
        End month of the seasonal window.
    output_csv : str, optional
        Destination path for the output CSV file. Default is
        ``"fire_recovery.csv"``.
    """
    fire_year = 2011
    lat, lon = 45.3311, -121.7331
    name = "Dollar Lake"
    veg_type = "conifer"
    start_year = fire_year - years_before
    end_year = fire_year + years_after

    print(f"Demo: {name} fire ({lat}, {lon}), year {fire_year}")
    print(f"Analysis window: {start_year}-{end_year}, months {month_start}-{month_end}\n")

    point = ee.Geometry.Point([lon, lat])
    buffer = point.buffer(500)

    collection = get_merged_collection(buffer, start_year, end_year,
                                       month_start, month_end)
    composites = annual_composites(collection, start_year, end_year,
                                   month_start, month_end)
    ts = extract_values_batch(composites, point)

    print(f"{'Year':>6}  {'NBR':>8}  {'NDVI':>8}")
    print("-" * 28)
    for r in ts:
        nbr = f"{r['NBR']:.4f}" if r["NBR"] is not None else "   N/A"
        ndvi = f"{r['NDVI']:.4f}" if r["NDVI"] is not None else "   N/A"
        print(f"{r['year']:>6}  {nbr:>8}  {ndvi:>8}")

    print()
    for index in ("NBR", "NDVI"):
        rec = compute_recovery(ts, fire_year, index)
        print(f"{index} Recovery:")
        print(f"  Pre-fire mean:  {rec['pre_fire']}")
        print(f"  Post-fire min:  {rec['post_fire_min']}")
        print(f"  Latest value:   {rec['latest']} ({rec['latest_year']})")
        print(f"  Drop (delta):   {rec['delta']}")
        print(f"  Recovery ratio: {rec['recovery_ratio']}")
        print()

    # Write CSV
    recovery_nbr = compute_recovery(ts, fire_year, "NBR")
    recovery_ndvi = compute_recovery(ts, fire_year, "NDVI")

    rows = []
    for r in ts:
        rows.append({
            "name": name,
            "veg_type": veg_type,
            "fire_year": fire_year,
            **r,
            **{f"NBR_{k}": v for k, v in recovery_nbr.items()},
            **{f"NDVI_{k}": v for k, v in recovery_ndvi.items()},
        })

    keys = rows[0].keys()
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Results written to {output_csv}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Parse command-line arguments and run the analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze vegetation recovery in fire-affected areas "
                    "using Landsat Collection 2 time series (NBR & NDVI)."
    )
    parser.add_argument(
        "--fire-shapefile", "-f",
        help="Path to fire areas shapefile. Must contain 'fire_year' column, "
             "optionally 'veg_type' and 'name'. If omitted, runs a demo.",
    )
    parser.add_argument(
        "--output", "-o", default="fire_recovery.csv",
        help="Output CSV path (default: fire_recovery.csv)",
    )
    parser.add_argument(
        "--years-before", type=int, default=3,
        help="Years of pre-fire data to include (default: 3)",
    )
    parser.add_argument(
        "--years-after", type=int, default=10,
        help="Years of post-fire data to include (default: 10)",
    )
    parser.add_argument(
        "--month-start", type=int, default=7,
        help="Start month for annual window (default: 7 = July)",
    )
    parser.add_argument(
        "--month-end", type=int, default=8,
        help="End month for annual window (default: 8 = August)",
    )
    parser.add_argument(
        "--project", "-p",
        help="GEE cloud project ID (required for newer ee API versions)",
    )

    args = parser.parse_args()

    try:
        if args.project:
            ee.Initialize(project=args.project)
        else:
            ee.Initialize()
    except ee.EEException:
        print("Earth Engine not authenticated. Run: earthengine authenticate")
        sys.exit(1)

    if args.fire_shapefile:
        process_shapefile(
            args.fire_shapefile,
            args.years_before,
            args.years_after,
            args.month_start,
            args.month_end,
            args.output,
        )
    else:
        print("No shapefile provided — running demo.\n")
        run_demo(args.years_before, args.years_after,
                 args.month_start, args.month_end, args.output)


if __name__ == "__main__":
    main()
