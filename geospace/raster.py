import os
import numpy as np
from osgeo import gdal
from pathlib import Path
from geospace.projection import read_srs
from geospace._const import WGS84, CREATION, TYPE_MAP
from geospace.utils import rep_file, ds_name, context_file


def convert_uint8(ds, in_no_data=None, out_no_data=255):
    ds, ras = ds_name(ds)
    frist_band = ds.GetRasterBand(1)

    if (frist_band.DataType != gdal.GDT_Byte or
            frist_band.ReadAsArray(0, 0, 1, 1).dtype != np.int8):
        return ras

    if frist_band.GetNoDataValue() is not None:
        in_no_data = frist_band.GetNoDataValue()
    if in_no_data is None:
        raise (ValueError("in_no_data must be initialed"))

    option = gdal.WarpOptions(multithread=True,
                              creationOptions=CREATION,
                              srcNodata=in_no_data,
                              dstNodata=out_no_data,
                              outputType=gdal.GDT_Byte)
    out_file = rep_file(os.path.dirname(ras), ras)
    gdal.Warp(out_file, ras, options=option)

    ds = None
    gdal.GetDriverByName('GTiff').Delete(ras)
    os.rename(out_file, ras)

    return ras


def resample(ds, out_path, **kwargs):
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    option = gdal.WarpOptions(multithread=True,
                              creationOptions=CREATION,
                              resampleAlg=resample_alg,
                              **kwargs)
    gdal.Warp(out_file, ds, options=option)

    return out_file


def mosaic(ras_paths, out_path, **kwargs):
    ds = ras_paths[0]
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    separate = kwargs.pop('separate', False)
    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    ds = gdal.BuildVRT('/vsimem/Mosaic.vrt', ras_paths, separate=separate)

    option = gdal.WarpOptions(multithread=True,
                              creationOptions=CREATION,
                              resampleAlg=resample_alg,
                              **kwargs)
    ds_out = gdal.Warp(out_file, ds, options=option)

    if separate:
        # each raster only have one band in the mosaic
        band_names = (Path(p).stem for p in ras_paths)
        [ds_out.GetRasterBand(i + 1).SetDescription(band_name)
         for i, band_name in enumerate(band_names)]

    return out_file


def project_raster(ds, out_path, **kwargs):
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    # input SpatialReference
    in_srs = kwargs.pop('srcSRS', None)
    inSpatialRef = read_srs([ds, in_srs])

    # output SpatialReference
    out_srs = kwargs.pop('dstSRS', WGS84)
    outSpatialRef = read_srs(out_srs)

    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    option = gdal.WarpOptions(creationOptions=CREATION,
                              resampleAlg=resample_alg,
                              srcSRS=inSpatialRef,
                              dstSRS=outSpatialRef,
                              multithread=True, **kwargs)
    gdal.Warp(out_file, ds, options=option)

    return out_file


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

    srs = kwargs.pop('dstSRS', WGS84)
    option = gdal.WarpOptions(multithread=True,
                              dstSRS=read_srs(srs),
                              creationOptions=CREATION,
                              **kwargs)
    gdal.Warp(out_file, ds, options=option)

    return out_file


def tif_copy_assign(out_file, ds_eg, array, srs=None, no_data=None):
    if os.path.exists(out_file):
        return out_file
    ds_eg = ds_name(ds_eg)[0]

    if array.ndim == 2:
        array = array.reshape([1, *array.shape])
    if array.ndim != 3:
        raise (Exception('array dims must be 2 or 3'))

    # set nodata value
    if no_data is None:
        if ds_eg.GetRasterBand(1).GetNoDataValue() is not None:
            no_data = ds_eg.GetRasterBand(1).GetNoDataValue()
        else:
            raise (Exception('nodata must be passed'))
    if isinstance(array, np.ma.core.MaskedArray):
        array.set_fill_value(no_data)
        array = array.filled()

    ds = gdal.GetDriverByName('GTiff').Create(
        out_file, array.shape[2], array.shape[1], array.shape[0],
        TYPE_MAP[array.dtype.name], CREATION)

    # fill with array
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(no_data)
    ds.WriteArray(array)

    # set geotransform
    trans = ds_eg.GetGeoTransform()
    ds.SetGeoTransform(tuple(trans))

    # set SpatialReference
    ds.SetSpatialRef(read_srs([srs, ds_eg]))

    return out_file
