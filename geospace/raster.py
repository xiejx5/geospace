import os
import numpy as np
from osgeo import gdal, osr, ogr
from geospace._const import CREATION, TYPE_MAP
from geospace.utils import block_write, rep_file, ds_name, context_file


def _map_no_data(in_data):
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
        block_write(ds, band, band_out, _map_no_data)

    ds = None
    ds_out = None
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
    gdal.BuildVRT('/vsimem/Mosaic.vrt', ras_paths, separate=separate)
    ds = gdal.Open('/vsimem/Mosaic.vrt')

    # set no data
    if ds.GetRasterBand(1).GetNoDataValue() is not None:
        no_data = ds.GetRasterBand(1).GetNoDataValue()
    else:
        no_data = kwargs.pop('srcNodata', None)

    option = gdal.WarpOptions(multithread=True,
                              dstNodata=no_data,
                              srcNodata=no_data,
                              creationOptions=CREATION,
                              resampleAlg=resample_alg,
                              **kwargs)
    gdal.Warp(out_file, ds, options=option)

    return out_file


def project_raster(ds, out_path, **kwargs):
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    # input SpatialReference
    in_srs = kwargs.pop('srcSRS', None)
    inSpatialRef = ds.GetSpatialRef()
    if not inSpatialRef:
        if in_srs is not None:
            if isinstance(in_srs, osr.SpatialReference):
                inSpatialRef = in_srs
            else:
                inSpatialRef = osr.SpatialReference()
                inSpatialRef.ImportFromProj4(in_srs)
        else:
            raise(ValueError("srcSRS must be set"))

    # output SpatialReference
    out_srs = kwargs.pop('dstSRS', "+proj=longlat +datum=WGS84 +ellps=WGS84")
    if isinstance(out_srs, osr.SpatialReference):
        outSpatialRef = out_srs
    elif os.path.isfile(out_srs):
        ext = os.path.splitext(os.path.basename(out_srs))[1]
        if ext == '.tif':
            ds_temp = gdal.Open(out_srs)
            outSpatialRef = ds_temp.GetSpatialRef()
            ds_temp = None
        else:
            ds_temp = ogr.Open(out_srs)
            outSpatialRef = ds_temp.GetLayer().GetSpatialRef()
            ds_temp = None
    else:
        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromProj4(out_srs)

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

    SpatialRef = osr.SpatialReference()
    SpatialRef.ImportFromProj4("+proj=longlat +datum=WGS84 +ellps=WGS84")
    srs = kwargs.pop('dstSRS', SpatialRef)
    option = gdal.WarpOptions(multithread=True,
                              dstSRS=srs,
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
        raise(Exception('array must be 2 dims or 3 dims'))

    # set nodata value
    if no_data is None:
        if ds_eg.GetRasterBand(1).GetNoDataValue() is not None:
            no_data = ds_eg.GetRasterBand(1).GetNoDataValue()
        else:
            raise(Exception('nodata must be passed'))
    if isinstance(array, np.ma.core.MaskedArray):
        array.set_fill_value(no_data)
        array = array.filled()

    ds = gdal.GetDriverByName('GTiff').Create(
        out_file, array.shape[2], array.shape[1], array.shape[0],
        TYPE_MAP[array.dtype.name], CREATION)

    # fill with array
    for i in range(1, 1 + ds.RasterCount):
        band = ds.GetRasterBand(i)
        band.SetNoDataValue(no_data)
        band.WriteArray(array[i - 1])

    # set geotransform
    trans = ds_eg.GetGeoTransform()
    ds.SetGeoTransform(tuple(trans))

    # set SpatialReference
    if srs is None:
        if ds_eg.GetProjection():
            inSpatialRef = osr.SpatialReference(wkt=ds_eg.GetProjection())
        else:
            raise(Exception('SpatialReference must be passed'))
    else:
        if isinstance(srs, osr.SpatialReference):
            inSpatialRef = srs
        else:
            inSpatialRef = osr.SpatialReference()
            inSpatialRef.ImportFromProj4(srs)
    ds.SetProjection(inSpatialRef.ExportToWkt())

    band = None
    ds = None

    return out_file
