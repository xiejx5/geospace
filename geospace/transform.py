import os
import json
import numpy as np
from osgeo import gdal

__all__ = ['context_file', 'ds_name']


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
