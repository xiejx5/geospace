import os
import gdal
from .transform import ds_name, context_file
from ._const import CREATION, CONFIG

__all__ = ['mosaic']


def mosaic(ras_paths, out_path, **kwargs):
    ds = ras_paths[0]
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    separate = kwargs.pop('separate', False)
    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    gdal.BuildVRT('/vsimem/Mosaic.vrt', ras_paths, separate=separate)
    ds = gdal.Open('/vsimem/Mosaic.vrt')
    option = gdal.WarpOptions(multithread=True, options=CONFIG,
                              creationOptions=CREATION, **kwargs,
                              resampleAlg=resample_alg)
    gdal.Warp(out_file, ds, options=option)

    return out_file
