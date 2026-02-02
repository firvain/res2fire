# res2fire

Fetch annual NBR and NDVI time series from Landsat Collection 2 for fire-affected areas.

Extracts annual median Normalized Burn Ratio (NBR) and Normalized Difference Vegetation Index (NDVI) values from Landsat 5, 7, 8, and 9 surface reflectance data via Google Earth Engine. Output is a CSV of raw time series values — analysis is performed separately.

## Prerequisites

- Python 3.10+
- [uv](https://docs.astral.sh/uv/) package manager
- A [Google Earth Engine](https://earthengine.google.com/) account with a cloud project

## Setup

```bash
# Clone and install dependencies
git clone <repo-url>
cd res2fire
uv sync

# Authenticate with Google Earth Engine (one-time)
uv run earthengine authenticate
```

## Usage

### Demo mode

Fetch a time series for the Dollar Lake fire location (Mt. Hood, Oregon):

```bash
uv run res2fire.py --project YOUR_GEE_PROJECT
```

### With a shapefile

```bash
uv run res2fire.py \
  --fire-shapefile fires.shp \
  --output results.csv \
  --project YOUR_GEE_PROJECT
```

### All options

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--fire-shapefile` | `-f` | *(demo mode)* | Path to fire areas vector file |
| `--output` | `-o` | `fire_timeseries.csv` | Output CSV path |
| `--start-year` | | `1985` | First year of the time series |
| `--end-year` | | `2024` | Last year of the time series |
| `--month-start` | | `6` | Start month of seasonal window (1-12) |
| `--month-end` | | `8` | End month of seasonal window (1-12) |
| `--project` | `-p` | | GEE cloud project ID |

### Example

Fetch time series from 2000 to 2023 using a July-September window:

```bash
uv run res2fire.py \
  -f fires.shp \
  -o timeseries.csv \
  --start-year 2000 \
  --end-year 2023 \
  --month-start 7 \
  --month-end 9 \
  -p my-gee-project
```

## Input shapefile format

Any OGR-readable vector format (Shapefile, GeoPackage, GeoJSON, etc.) in **WGS84 (EPSG:4326)**.

### Required

| Column | Type | Description |
|--------|------|-------------|
| geometry | Point, Polygon, or MultiPolygon | Location or perimeter of the area of interest |

### Optional columns

| Column | Type | Default | Description |
|--------|------|---------|-------------|
| `veg_type` | str | `unknown` | Vegetation type: `broadleaf`, `conifer`, or `macchia` |
| `name` or `id` | str | row index | Human-readable identifier |

### Notes

- **Points** are buffered by 500 m for image filtering; the index value is sampled at the exact point location.
- **Polygons / MultiPolygons** are used directly — the script computes the spatial mean of index values over the entire area.
- If the shapefile from partners uses different column names, rename them to match before running.

### Example attribute table

| name | veg_type | geometry |
|------|----------|----------|
| Fire_A | broadleaf | POINT(23.5, 39.2) |
| Fire_B | conifer | POLYGON(...) |
| Fire_C | macchia | POLYGON(...) |

## Output

The output CSV contains one row per site per year:

| Column | Description |
|--------|-------------|
| `name` | Site identifier |
| `veg_type` | Vegetation type |
| `year` | Observation year |
| `NBR` | Normalized Burn Ratio (annual median, June-August) |
| `NDVI` | Normalized Difference Vegetation Index (annual median, June-August) |

## How it works

1. **Data loading** — Merges Landsat 5 TM, 7 ETM+, 8 OLI, and 9 OLI-2 Collection 2 Level-2 surface reflectance imagery from GEE.
2. **Preprocessing** — Applies USGS scale factors, masks clouds/shadows via QA_PIXEL, and renames bands to a common naming scheme.
3. **Index computation** — Calculates NBR and NDVI per image.
4. **Annual composites** — Produces a pixel-wise median for the seasonal window (default June-August) of each year.
5. **Value extraction** — Computes spatial mean of indices over each geometry in a single batched GEE request.

## References

- USGS, [Landsat Collection 2 Level-2 Science Products](https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products)
- Key, C.H. and Benson, N.C. (2006), "Landscape Assessment (LA)", FIREMON: Fire Effects Monitoring and Inventory System, RMRS-GTR-164-CD
