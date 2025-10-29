import os
import numpy as np
from osgeo import gdal
from collections.abc import Iterable
from geospace._const import WGS84, CREATION
from geospace.projection import read_srs, coord_trans


def block_write(datasets, map_fun, *args, **kwargs):
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
    ds = ds_name(ds)[0]
    trans = ds.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    col, row = np.linalg.solve(a, b) - 0.5
    return np.round(col).astype(int), np.round(row).astype(int)


def imagexy2geo(ds, row, col):
    """
    row, col to centroid lon, lat
    """
    ds = ds_name(ds)[0]
    trans = ds.GetGeoTransform()
    px = trans[0] + (col + 0.5) * trans[1] + (row + 0.5) * trans[2]
    py = trans[3] + (col + 0.5) * trans[4] + (row + 0.5) * trans[5]
    return px, py


def meshgrid(ds, geo_srs=WGS84):
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
    if isinstance(ds, gdal.Dataset):
        ras = ds.GetDescription()
    else:
        ras = str(ds)
        ds = gdal.Open(ras)

    return ds, ras


def get_extent(layer):
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
