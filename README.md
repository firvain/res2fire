# res2fire

Analyze vegetation recovery in fire-affected areas using Landsat Collection 2 time series.

Computes annual NBR (Normalized Burn Ratio) and NDVI time series from Landsat 5, 7, 8, and 9 surface reflectance data via Google Earth Engine, then quantifies post-fire recovery capacity.

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

Run the built-in demo on the 2011 Dollar Lake fire (Mt. Hood, Oregon):

```bash
uv run res2fire.py --project YOUR_GEE_PROJECT
```

### With a fire shapefile

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
| `--output` | `-o` | `fire_recovery.csv` | Output CSV path |
| `--years-before` | | `3` | Pre-fire years to include |
| `--years-after` | | `10` | Post-fire years to include |
| `--month-start` | | `7` | Start month of seasonal window (1-12) |
| `--month-end` | | `8` | End month of seasonal window (1-12) |
| `--project` | `-p` | | GEE cloud project ID |

### Example

Analyze 5 years before and 15 years after each fire, using a June-September window:

```bash
uv run res2fire.py \
  -f fires.shp \
  -o recovery_long.csv \
  --years-before 5 \
  --years-after 15 \
  --month-start 6 \
  --month-end 9 \
  -p my-gee-project
```

## Input shapefile format

Any OGR-readable vector format (Shapefile, GeoPackage, GeoJSON, etc.) in **WGS84 (EPSG:4326)**.

### Required columns

| Column | Type | Description |
|--------|------|-------------|
| `fire_year` | int | Year the fire occurred |
| geometry | Point, Polygon, or MultiPolygon | Location or perimeter of the burned area |

### Optional columns

| Column | Type | Default | Description |
|--------|------|---------|-------------|
| `veg_type` | str | `unknown` | Vegetation type: `broadleaf`, `conifer`, or `macchia` |
| `name` or `id` | str | row index | Human-readable fire identifier |

### Notes

- **Points** are buffered by 500 m for image filtering; the index value is sampled at the exact point location.
- **Polygons / MultiPolygons** are used directly — the script computes the spatial mean of index values over the entire area.
- If the shapefile from partners uses different column names, rename them to match before running.

### Example attribute table

| name | fire_year | veg_type | geometry |
|------|-----------|----------|----------|
| Fire_A | 2017 | broadleaf | POINT(23.5, 39.2) |
| Fire_B | 2015 | conifer | POLYGON(...) |
| Fire_C | 2019 | macchia | POLYGON(...) |

## Output

The output CSV contains one row per fire site per year, with the following columns:

| Column | Description |
|--------|-------------|
| `name` | Fire identifier |
| `veg_type` | Vegetation type |
| `fire_year` | Year of fire |
| `year` | Observation year |
| `NBR` | Normalized Burn Ratio for that year |
| `NDVI` | Normalized Difference Vegetation Index for that year |
| `NBR_pre_fire` | Mean NBR over 3 pre-fire years |
| `NBR_post_fire_min` | Minimum post-fire NBR |
| `NBR_latest` | Most recent valid NBR value |
| `NBR_latest_year` | Year of the latest NBR value |
| `NBR_delta` | Fire-induced NBR drop |
| `NBR_recovery_ratio` | Recovery ratio (1.0 = full recovery) |
| `NDVI_*` | Same metrics for NDVI |

## How it works

1. **Data loading** -- Merges Landsat 5 TM, 7 ETM+, 8 OLI, and 9 OLI-2 Collection 2 Level-2 surface reflectance imagery from GEE.
2. **Preprocessing** -- Applies USGS scale factors, masks clouds/shadows via QA_PIXEL, and renames bands to a common naming scheme.
3. **Index computation** -- Calculates NBR and NDVI per image.
4. **Annual composites** -- Produces a pixel-wise median for the seasonal window (default July-August) of each year.
5. **Value extraction** -- Computes spatial mean of indices over each fire geometry in a single batched GEE request.
6. **Recovery metrics** -- Compares pre-fire baseline (3-year mean) with post-fire trajectory to derive a recovery ratio.

### Recovery ratio

```
RR = (V_latest - V_post_min) / (V_pre - V_post_min)
```

- `RR = 1.0` -- vegetation has fully recovered to pre-fire levels
- `RR > 1.0` -- vegetation exceeds pre-fire levels
- `RR < 1.0` -- recovery is still in progress

## References

- USGS, [Landsat Collection 2 Level-2 Science Products](https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products)
- Key, C.H. and Benson, N.C. (2006), "Landscape Assessment (LA)", FIREMON: Fire Effects Monitoring and Inventory System, RMRS-GTR-164-CD
