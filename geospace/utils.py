import os
import json
import numpy as np
from osgeo import gdal
from collections.abc import Iterable
from geospace._const import CREATION
from geospace.projection import read_srs, coord_trans


def block_write(in_ds, in_bands, out_band, map_fun,
                map_args=None, map_kwargs=None):
    if isinstance(in_bands, Iterable):
        _block_write_multi(in_ds, in_bands, out_band, map_fun,
                           map_args=map_args, map_kwargs=map_kwargs)
    else:
        _block_write_one(in_ds, in_bands, out_band, map_fun,
                         map_args=map_args, map_kwargs=map_kwargs)


def _block_write_one(in_ds, in_band, out_band, map_fun,
                     map_args=None, map_kwargs=None):
    block_xsize, block_ysize = in_band.GetBlockSize()
    for b_y, yoff in enumerate(range(0, in_ds.RasterYSize, block_ysize)):
        for b_x, xoff in enumerate(range(0, in_ds.RasterXSize, block_xsize)):
            win_xsize, win_ysize = in_band.GetActualBlockSize(b_x, b_y)

            in_data = in_band.ReadAsArray(
                xoff=xoff, yoff=yoff, win_xsize=win_xsize, win_ysize=win_ysize)
            in_data = np.ma.masked_array(in_data,
                                         mask=in_data == in_band.GetNoDataValue(),
                                         fill_value=out_band.GetNoDataValue())

            if map_args and map_kwargs:
                out_data = map_fun(in_data, *map_args, **map_kwargs)
            elif map_args:
                out_data = map_fun(in_data, *map_args)
            elif map_kwargs:
                out_data = map_fun(in_data, **map_kwargs)
            else:
                out_data = map_fun(in_data)

            out_band.WriteArray(out_data, xoff=xoff, yoff=yoff)


def _block_write_multi(in_ds, in_bands, out_band, map_fun,
                       map_args=None, map_kwargs=None):
    first_band = in_bands[0]
    block_xsize, block_ysize = first_band.GetBlockSize()
    for b_y, yoff in enumerate(range(0, in_ds.RasterYSize, block_ysize)):
        for b_x, xoff in enumerate(range(0, in_ds.RasterXSize, block_xsize)):
            win_xsize, win_ysize = first_band.GetActualBlockSize(b_x, b_y)
            no_data = out_band.GetNoDataValue()

            in_datas = []
            for in_band in in_bands:
                in_data = in_band.ReadAsArray(
                    xoff=xoff, yoff=yoff, win_xsize=win_xsize, win_ysize=win_ysize)
                in_data = np.ma.masked_array(in_data,
                                             mask=in_data == in_band.GetNoDataValue(),
                                             fill_value=no_data)
                in_datas.append(in_data)

            if map_args and map_kwargs:
                out_data = map_fun(in_datas, *map_args, **map_kwargs)
            elif map_args:
                out_data = map_fun(in_datas, *map_args)
            elif map_kwargs:
                out_data = map_fun(in_datas, **map_kwargs)
            else:
                out_data = map_fun(in_datas)

            out_band.WriteArray(out_data, xoff=xoff, yoff=yoff)


def rep_file(cache_dir, filename):
    prefix, extension = os.path.splitext(os.path.basename(filename))
    file_path = os.path.join(cache_dir, prefix + extension)
    if os.path.exists(file_path):
        i = 1
        while True:
            file_path = os.path.join(
                cache_dir, prefix + '(' + str(i) + ')' + extension)
            if os.path.exists(file_path):
                i += 1
            else:
                break
    return file_path


def geo2imagexy(ds, x, y):
    trans = ds.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    col, row = np.linalg.solve(a, b) - 0.5
    return int(round(col)), int(round(row))


def imagexy2geo(dataset, row, col):
    '''
    row, col to centroid lon, lat
    '''
    trans = dataset.GetGeoTransform()
    px = trans[0] + (col + 0.5) * trans[1] + (row + 0.5) * trans[2]
    py = trans[3] + (col + 0.5) * trans[4] + (row + 0.5) * trans[5]
    return px, py


def meshgrid(ds, geo_srs="+proj=longlat +datum=WGS84 +ellps=WGS84"):
    ds, _ = ds_name(ds)
    col, row = np.meshgrid(np.arange(ds.RasterXSize), np.arange(ds.RasterYSize))
    x, y = imagexy2geo(ds, row, col)

    if 'PROJCS' not in ds.GetProjection():
        return x, y

    trans = coord_trans(ds, geo_srs)
    lonlat = np.array(trans.TransformPoints(np.concatenate([x.reshape(-1, 1), y.reshape(-1, 1)], axis=1)))
    lon = lonlat[:, 0].reshape(x.shape)
    lat = lonlat[:, 1].reshape(y.shape)
    return lon, lat


def context_file(ras, out_path):
    ext = os.path.splitext(os.path.basename(out_path))[1]
    if ext != '.tif':
        if not os.path.isdir(out_path):
            os.makedirs(out_path)
        out_file = os.path.join(out_path, os.path.splitext(
            os.path.basename(ras))[0] + '.tif')
    else:
        out_file = out_path

    return out_file


def ds_name(ds):
    if isinstance(ds, str):
        ras = ds
        ds = gdal.Open(ras)
    else:
        ras = ds.GetDescription()

    return ds, ras


def get_extent(layer):
    r = layer.GetSpatialRef()
    crs = json.loads(r.ExportToPROJJSON())
    axis = [crs['coordinate_system']['axis'][i]['direction'].upper()
            for i in range(r.GetAxesCount())]
    order = ['EAST', 'WEST', 'NORTH', 'SOUTH']
    order_dict = dict(zip(order, range(4)))
    axis_code = np.array([order_dict[i] for i in axis])
    axis_sort = axis_code.argsort()
    data_map = [abs(i) - 1 for i in r.GetDataAxisToSRSAxisMapping()]
    extent = np.reshape(layer.GetExtent(), [-1, 2])[data_map]
    xy_extent = np.reshape(extent[axis_sort], [-1])
    return tuple(xy_extent)


def zeros_tif(out_file, x_size, y_size, n_band,
              data_type, trans, srs, no_data=2):

    ds = gdal.GetDriverByName('GTiff').Create(
        out_file, x_size, y_size, n_band, data_type, CREATION)

    # fill with 0 and set no data
    band = ds.GetRasterBand(1)
    band.Fill(0)
    band.SetNoDataValue(no_data)

    # set geotransform
    ds.SetGeoTransform(tuple(trans))

    # set SpatialReference
    ds.SetSpatialRef(read_srs(srs))

    return out_file
