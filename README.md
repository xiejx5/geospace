# GeoSpace

**A Python library for geospatial processing based on GDAL.**

GeoSpace provides practical tools for raster and vector processing, including zonal statistics, reprojection, resampling, mosaicking, clipping, rasterization, and Google Earth Engine workflows.

## Installation

```bash
uv add geospace
```

GeoSpace wheels include the GDAL runtime and Python bindings, so no separate GDAL installation is required.

Both imports work directly:

```python
import geospace
from osgeo import gdal
```

The bundled GDAL wheels are built and tested separately in [gdal-wheels](https://github.com/xiejx5/gdal-wheels).

> Do not install a separate `gdal` package in the same environment.

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
- GRIB to GeoTIFF conversion
- Coordinate and spatial calculations
- Optional Google Earth Engine integration
