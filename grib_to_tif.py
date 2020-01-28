import os
from osgeo import gdal, osr
from .transform import ds_name, context_file
from ._const import CREATION, CONFIG


__all__ = ['grib_to_tif']


def grib_to_tif(ds, out_path=None, **kwargs):
    ds, ras = ds_name(ds)

    if os.path.splitext(os.path.basename(ras))[1] != '.grib':
        return

    if out_path:
        out_file = context_file(ras, out_path)
    else:
        out_file = os.path.join(os.path.dirname(ras), os.path.splitext(
            os.path.basename(ras))[0] + '.tif')

    if os.path.exists(out_file):
        return out_file

    proj = "+proj=longlat +datum=WGS84 +ellps=WGS84"
    SpatialRef = osr.SpatialReference()
    SpatialRef.ImportFromProj4(proj)
    srs = kwargs.pop('dstSRS', SpatialRef)
    option = gdal.WarpOptions(multithread=True, options=CONFIG,
                              creationOptions=CREATION, **kwargs,
                              dstSRS=srs)
    gdal.Warp(out_file, ds, options=option)

    return out_file
