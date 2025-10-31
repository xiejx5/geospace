# GeoSpace

**A Python library for geospatial processing based on GDAL.**

GeoSpace is a powerful and efficient Python library designed to simplify geospatial data processing. It provides a comprehensive set of tools for working with raster and vector data, leveraging the capabilities of GDAL.

## Installation

You can install GeoSpace using pip:

```bash
pip install geospace -U
```

## Features

### High-Performance Zonal Statistics
GeoSpace excels at performing zonal statistics, offering highly optimized functions for calculating statistics of raster values within polygons.

- **Area-Weighted Averaging:** The `basin_average` function calculates the area-weighted average of raster values for each polygon in a shapefile. It is optimized for memory efficiency and parallel processing, making it suitable for large datasets and high-performance computing environments, including Slurm clusters.

#### `basin_average` Example
Here's an example of how to use the `basin_average` function to calculate the average precipitation and temperature for a set of basins:

```python
import geospace as gs

# Paths to your rasters and shapefile
rasters = ['path/to/P.tif', 'path/to/T.tif']
basins_shp = 'path/to/basins.shp'

# Calculate the average values for each basin
df_stats = gs.basin_average(rasters, basins_shp)

# Rasters as rows and basins as columns
print(df_stats)
```

### Raster Processing
- **Reprojection:** Easily reproject rasters to different coordinate systems.
- **Resampling:** Resample rasters to new resolutions using various algorithms.
- **Mosaicking:** Combine multiple rasters into a single mosaic.
- **Clipping:** Clip rasters to the extent of a shapefile.
- **Nodata Filling:** Fill nodata values in rasters using various interpolation methods.
- **Data Type Conversion:** Convert rasters to different data types (e.g., UInt8).
- **GRIB to GeoTIFF:** Convert GRIB files to GeoTIFF format.

### Vector Processing
- **Buffering:** Create buffers around vector features.
- **Reprojection:** Reproject shapefiles to different coordinate systems.
- **Filtering:** Filter shapefiles based on attribute queries.
- **Polygonization:** Convert rasters to polygon shapefiles.
- **Rasterization:** Convert shapefiles to raster datasets.

### Google Earth Engine Integration
- **Data Export:** Export Google Earth Engine images to GeoTIFF or CSV format.
- **SoilGrids:** Download and preprocess SoilGrids data.
- **Wind Data:** Calculate wind speed from u and v components.
- **Time Series Analysis:** Group image collections by month and calculate seasonality indices.

### Utilities
- **Coordinate Conversion:** Convert between geographic and image coordinates.
- **Spatial Calculations:** Calculate grid cell areas and distances.
- **File Handling:** Manage file paths and avoid name collisions.
