import os
from osgeo import gdal
from .transform import ds_name, context_file
from ._const import CREATION, CONFIG


__all__ = ['resample']


def resample(ds, out_path, **kwargs):
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    option = gdal.WarpOptions(multithread=True, options=CONFIG,
                              creationOptions=CREATION, **kwargs,
                              resampleAlg=resample_alg)
    gdal.Warp(out_file, ds, options=option)

    return out_file
