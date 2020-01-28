import os
import numpy as np
from osgeo import gdal
from ._const import CREATION
from .block_write import block_write
from .transform import rep_file, ds_name

__all__ = ['convert_uint8']


def map_no_data(in_data):
    return in_data.filled()


def convert_uint8(ds, in_no_data=None, out_no_data=255):
    ds, ras = ds_name(ds)
    frist_band = ds.GetRasterBand(1)

    if (frist_band.DataType != gdal.GDT_Byte or
            frist_band.ReadAsArray(0, 0, 1, 1).dtype != np.int8):
        return ras

    if frist_band.GetNoDataValue() is not None:
        in_no_data = frist_band.GetNoDataValue()
    if in_no_data is None:
        raise(ValueError("in_no_data must be initialed"))

    option = gdal.WarpOptions(multithread=True,
                              creationOptions=CREATION,
                              outputType=gdal.GDT_Byte)
    out_file = rep_file(os.path.dirname(ras), ras)
    ds_out = gdal.Warp(out_file, ras, options=option)

    for i in range(1, 1 + ds.RasterCount):
        band = ds.GetRasterBand(i)
        band_out = ds_out.GetRasterBand(i)
        if band.GetNoDataValue() is None:
            band.SetNoDataValue(in_no_data)
        band_out.SetNoDataValue(out_no_data)
        block_write(ds, band, band_out, map_no_data)

    ds = None
    ds_out = None
    gdal.GetDriverByName('GTiff').Delete(ras)
    os.rename(out_file, ras)

    return ras
