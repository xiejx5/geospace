# GeoSpace

A lightweight Python library for geospatial processing based on GDAL.

GeoSpace provides practical tools for raster and vector processing, including zonal statistics, reprojection, resampling, mosaicking, clipping, rasterization, and optional Google Earth Engine workflows.

## Installation

```bash
uv add geospace
```

GeoSpace uses [`gdal-wheel`](https://github.com/xiejx5/gdal-wheel) as its GDAL dependency, so GDAL does not need to be compiled as part of GeoSpace.

```python
import geospace
from osgeo import gdal
```

## Example

```python
import geospace as gs

rasters = ["precipitation.tif", "temperature.tif"]
basins = "basins.shp"

stats = gs.reduce(rasters, basins)
print(stats)
```

## Features

- High-performance zonal statistics
- Raster reprojection, resampling, mosaicking, clipping, and nodata filling
- Raster and vector conversion
- Shapefile projection, filtering, buffering, and rasterization
- Coordinate and spatial calculations
- Optional Google Earth Engine integration

## Optional dependencies

Google Earth Engine support:

```bash
uv add "geospace[gee]"
```

Mask-related features:

```bash
uv add "geospace[mask]"
```

## Development

```bash
git clone https://github.com/xiejx5/geospace.git
cd geospace
uv sync
```

GeoSpace is a pure-Python package. GDAL wheels are built and published separately by [`gdal-wheel`](https://github.com/xiejx5/gdal-wheel).
