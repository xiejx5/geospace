import os
import numpy as np
from osgeo import gdal
from collections.abc import Iterable
from geospace._const import WGS84, CREATION
from geospace.projection import read_srs, coord_trans


def block_write(datasets, map_fun, *args, **kwargs):
    """Applies a function to raster datasets block by block and writes the result.

    This function is designed to process large raster datasets that do not fit into memory.

    Args:
        datasets (list or gdal.Dataset): A list of GDAL datasets or a single dataset.
        map_fun (function): The function to apply to each block.
        *args: Positional arguments to pass to the map_fun.
        **kwargs: Keyword arguments to pass to the map_fun.
    """
    import psutil

    if not isinstance(datasets, Iterable):
        datasets = [datasets]
    datasets = [ds_name(ds)[0] for ds in datasets]
    ds = datasets[0]
    n_x, n_y, band = ds.RasterXSize, ds.RasterYSize, ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    ratio = int(
        n_y
        * n_x
        * ds.RasterCount
        * len(datasets)
        / (
            psutil.virtual_memory().available
            * 0.5
            / (2 + ds.ReadAsArray(0, 0, 1, 1).dtype.itemsize)
        )
        + 1
    )
    if ratio == 1:
        arrs = []
        for d in datasets:
            arr = d.ReadAsArray()
            mask = (arr == d.GetRasterBand(1).GetNoDataValue(),)
            arrs.append(np.ma.masked_array(arr, mask=mask, fill_value=nodata))
        arr = map_fun(arrs, *args, **kwargs)
        ds.WriteArray(arr)
        del arr, arrs, mask
    else:
        b_xsize, b_ysize = band.GetBlockSize()
        xsize = max(n_x // ratio + 1, b_xsize)
        ysize = max(n_y // ratio + 1, b_ysize)
        xoffs, yoffs = range(0, n_x, xsize), range(0, n_y, ysize)
        for xoff in xoffs:
            for yoff in yoffs:
                win_xsize = min(xsize, n_x - xoff)
                win_ysize = min(ysize, n_y - yoff)

                arrs = []
                for d in datasets:
                    arr = d.ReadAsArray(
                        xoff=xoff, yoff=yoff, xsize=win_xsize, ysize=win_ysize
                    )
                    mask = (arr == d.GetRasterBand(1).GetNoDataValue(),)
                    arrs.append(np.ma.masked_array(arr, mask=mask, fill_value=nodata))
                part_args = (
                    a[yoff : yoff + win_ysize, xoff : xoff + win_xsize] for a in args
                )
                part_kwargs = {
                    k: v[yoff : yoff + win_ysize, xoff : xoff + win_xsize]
                    for k, v in kwargs.items()
                }
                arr = map_fun(arrs, *part_args, **part_kwargs)
                ds.WriteArray(arr, xoff=xoff, yoff=yoff)
                del arr, arrs, mask, part_args, part_kwargs


def rep_file(cache_dir, filename):
    """Generates a new file path in a directory, avoiding name collisions.

    If a file with the same name already exists, a number is appended to the
    base name (e.g., 'file(1).txt').

    Args:
        cache_dir (str): The directory to save the file in.
        filename (str): The original filename.

    Returns:
        str: The new, non-colliding file path.
    """
    prefix, extension = os.path.splitext(os.path.basename(filename))
    file_path = os.path.join(cache_dir, prefix + extension)
    if os.path.exists(file_path):
        i = 1
        while True:
            file_path = os.path.join(cache_dir, prefix + '(' + str(i) + ')' + extension)
            if os.path.exists(file_path):
                i += 1
            else:
                break
    return file_path


def rep_name(rasters, sort_idxs=None):
    """Generates descriptive names for raster bands.

    Args:
        rasters (list): A list of raster file paths.
        sort_idxs (np.ndarray, optional): An array of indices to sort the rasters by.
                                          Defaults to None.

    Returns:
        tuple: A tuple containing:
            - names (np.ndarray): An array of band names.
            - s (np.ndarray): The start index of each raster's bands.
            - t (np.ndarray): The end index of each raster's bands.
            - (optional) np.ndarray: The sorted indices.
    """
    n_bands = [gdal.Open(ras).RasterCount for ras in rasters]
    t = np.cumsum(n_bands)
    s = np.roll(t, 1)
    s[0] = 0

    names = np.zeros(t[-1], dtype='object')
    for i, ras in enumerate(rasters):
        string = os.path.splitext(os.path.basename(ras))[0]
        names[s[i] : t[i]] = np.strings.add(
            string, np.strings.mod('%d', np.arange(n_bands[i]))
        )
        names[s[i]] = string

    if sort_idxs is None:
        return names, s, t

    relative_loc = np.zeros(t[-1], dtype=int)
    idxs = sort_idxs * np.power(10, int(np.ceil(np.log10(np.max(n_bands)))))
    for i, _ in enumerate(rasters):
        relative_loc[s[i] : t[i]] = idxs[i] + np.arange(n_bands[i])
    return names, s, t, np.argsort(relative_loc)


def geo2imagexy(ds, x, y):
    """Converts geographic coordinates to image coordinates.

    Args:
        ds (gdal.Dataset or str): The raster dataset or its path.
        x (float): The geographic x-coordinate (longitude).
        y (float): The geographic y-coordinate (latitude).

    Returns:
        tuple: A tuple containing the column and row (col, row).
    """
    ds = ds_name(ds)[0]
    trans = ds.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    col, row = np.linalg.solve(a, b) - 0.5
    return np.round(col).astype(int), np.round(row).astype(int)


def imagexy2geo(ds, row, col):
    """Converts image coordinates to geographic coordinates of the cell's centroid.

    Args:
        ds (gdal.Dataset or str): The raster dataset or its path.
        row (int): The row index.
        col (int): The column index.

    Returns:
        tuple: A tuple containing the longitude and latitude (lon, lat).
    """
    ds = ds_name(ds)[0]
    trans = ds.GetGeoTransform()
    px = trans[0] + (col + 0.5) * trans[1] + (row + 0.5) * trans[2]
    py = trans[3] + (col + 0.5) * trans[4] + (row + 0.5) * trans[5]
    return px, py


def meshgrid(ds, geo_srs=WGS84):
    """Creates a meshgrid of geographic coordinates from a raster dataset.

    Args:
        ds (gdal.Dataset or str): The raster dataset or its path.
        geo_srs (str, optional): The output geographic spatial reference system.
                                 Defaults to WGS84.

    Returns:
        tuple: A tuple containing the longitude and latitude meshgrids (lon, lat).
    """
    ds, _ = ds_name(ds)
    col, row = np.meshgrid(np.arange(ds.RasterXSize), np.arange(ds.RasterYSize))
    x, y = imagexy2geo(ds, row, col)

    if 'PROJCS' not in ds.GetProjection():
        return x, y

    trans = coord_trans(ds, geo_srs)
    lonlat = np.array(
        trans.TransformPoints(
            np.concatenate([x.reshape(-1, 1), y.reshape(-1, 1)], axis=1)
        )
    )
    lon = lonlat[:, 0].reshape(x.shape)
    lat = lonlat[:, 1].reshape(y.shape)
    return lon, lat


def context_file(ras, out_path):
    """Generates an output file path based on a context raster.

    If out_path is a directory, the output filename is derived from the
    context raster's filename.

    Args:
        ras (str): The path to the context raster file.
        out_path (str): The output path (can be a file or a directory).

    Returns:
        str: The generated output file path.
    """
    out_path = str(out_path)
    ext = os.path.splitext(os.path.basename(out_path))[1]
    if ext != '.tif':
        if (not os.path.isdir(out_path)) and ('/vsimem' not in out_path):
            os.makedirs(out_path)
        out_file = os.path.join(
            out_path, os.path.splitext(os.path.basename(ras))[0] + '.tif'
        )
    else:
        out_file = out_path

    return out_file


def ds_name(ds):
    """Gets the GDAL dataset object and its name.

    Args:
        ds (gdal.Dataset or str): The GDAL dataset or its path.

    Returns:
        tuple: A tuple containing the GDAL dataset and its name (ds, name).
    """
    if isinstance(ds, gdal.Dataset):
        ras = ds.GetDescription()
    else:
        ras = str(ds)
        ds = gdal.Open(ras)

    return ds, ras


def get_extent(layer):
    """Gets the extent of an OGR layer.

    Args:
        layer (ogr.Layer): The OGR layer.

    Returns:
        tuple: The extent of the layer (x_min, y_min, x_max, y_max).
    """
    import json

    r = layer.GetSpatialRef()
    crs = json.loads(r.ExportToPROJJSON())
    axis = [
        crs['coordinate_system']['axis'][i]['direction'].upper()
        for i in range(r.GetAxesCount())
    ]
    order = ['EAST', 'WEST', 'NORTH', 'SOUTH']
    order_dict = dict(zip(order, range(4)))
    axis_code = np.array([order_dict[i] for i in axis])
    axis_sort = axis_code.argsort()
    data_map = [abs(i) - 1 for i in r.GetDataAxisToSRSAxisMapping()]
    extent = np.reshape(layer.GetExtent(), [-1, 2])[data_map]
    xy_extent = np.reshape(extent[axis_sort], [-1])
    return tuple(xy_extent)


def zeros_tif(out_file, x_size, y_size, n_band, data_type, trans, srs, nodata=2):
    """Creates a new GeoTIFF file filled with zeros.

    Args:
        out_file (str): The path to the output GeoTIFF file.
        x_size (int): The width of the raster in pixels.
        y_size (int): The height of the raster in pixels.
        n_band (int): The number of bands.
        data_type (int): The data type of the raster (e.g., gdal.GDT_Byte).
        trans (tuple): The geotransform.
        srs (osr.SpatialReference or str): The spatial reference system.
        nodata (int, optional): The nodata value. Defaults to 2.

    Returns:
        str: The path to the output GeoTIFF file.
    """
    ds = gdal.GetDriverByName('GTiff').Create(
        out_file, x_size, y_size, n_band, data_type, CREATION
    )

    # fill with 0 and set no data
    band = ds.GetRasterBand(1)
    band.Fill(0)
    band.SetNoDataValue(nodata)

    # set geotransform
    ds.SetGeoTransform(tuple(trans))

    # set SpatialReference
    ds.SetSpatialRef(read_srs(srs))

    return out_file
