import os
import numpy as np
from osgeo import gdal
from pathlib import Path
from geospace.projection import read_srs
from geospace._const import WGS84, CREATION, TYPE_MAP
from geospace.utils import rep_file, ds_name, context_file


def convert_uint8(ds, in_nodata=None, out_nodata=255):
    """Converts a raster dataset to UInt8 data type.

    Args:
        ds (gdal.Dataset or str): The input raster dataset or its path.
        in_nodata (float, optional): The input nodata value. If None, it is read from the dataset.
                                     Defaults to None.
        out_nodata (int, optional): The output nodata value. Defaults to 255.

    Returns:
        str: The path to the converted raster file.
    """
    ds, ras = ds_name(ds)
    frist_band = ds.GetRasterBand(1)

    if (
        frist_band.DataType != gdal.GDT_Byte
        or frist_band.ReadAsArray(0, 0, 1, 1).dtype != np.int8
    ):
        return ras

    if frist_band.GetNoDataValue() is not None:
        in_nodata = frist_band.GetNoDataValue()
    if in_nodata is None:
        raise (ValueError('in_nodata must be initialed'))

    option = gdal.WarpOptions(
        multithread=True,
        creationOptions=CREATION,
        srcNodata=in_nodata,
        dstNodata=out_nodata,
        outputType=gdal.GDT_Byte,
        warpOptions=['UNIFIED_SRC_NODATA=PARTIAL'],
    )
    out_file = rep_file(os.path.dirname(ras), ras)
    gdal.Warp(out_file, ras, options=option)

    ds = None
    gdal.GetDriverByName('GTiff').Delete(ras)
    os.rename(out_file, ras)

    return ras


def resample(ds, out_path, **kwargs):
    """Resamples a raster dataset to a new resolution.

    Args:
        ds (gdal.Dataset or str): The input raster dataset or its path.
        out_path (str): The path to the output resampled raster file.
        **kwargs: Additional arguments for gdal.WarpOptions().

    Returns:
        str: The path to the output resampled raster file.
    """
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    warp_options = kwargs.pop('warpOptions', []) + ['UNIFIED_SRC_NODATA=PARTIAL']
    option = gdal.WarpOptions(
        multithread=True,
        creationOptions=CREATION,
        resampleAlg=resample_alg,
        warpOptions=warp_options,
        **kwargs,
    )
    gdal.Warp(out_file, ds, options=option)

    return out_file


def mosaic(rasters, out_path, **kwargs):
    """Mosaics multiple raster datasets into a single file.

    Args:
        rasters (list or str): A list of raster file paths or a single path.
        out_path (str): The path to the output mosaic raster file.
        **kwargs: Additional arguments for gdal.WarpOptions().

    Returns:
        str: The path to the output mosaic raster file.
    """
    if isinstance(rasters, (str, Path)):
        rasters = [str(rasters)]
    else:
        rasters = [str(ras) for ras in rasters]

    ds, ras = ds_name(rasters[0])
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    # convert longitude [0, 360] to [-180, 180] when mosaic nc files
    all_nc = all(Path(f).suffix == '.nc' for f in rasters)
    separate = kwargs.pop('separate', True if all_nc else False)
    srcSRS = read_srs([ds, kwargs.pop('srcSRS', WGS84)])
    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    outputBounds = kwargs.pop('outputBounds', [-180, -90, 180, 90] if all_nc else None)
    warp_options = kwargs.pop('warpOptions', []) + ['UNIFIED_SRC_NODATA=PARTIAL']
    option = gdal.WarpOptions(
        multithread=True,
        srcSRS=srcSRS,
        creationOptions=CREATION,
        resampleAlg=resample_alg,
        outputBounds=outputBounds,
        warpOptions=warp_options,
        **kwargs,
    )

    ds = gdal.BuildVRT('/vsimem/Mosaic.vrt', rasters, separate=separate)
    ds_out = gdal.Warp(out_file, ds, options=option)

    if separate:
        # each raster only have one band in the mosaic
        band_names = (Path(p).stem for p in rasters)
        [
            ds_out.GetRasterBand(i + 1).SetDescription(band_name)
            for i, band_name in enumerate(band_names)
        ]

    return out_file


def project_raster(ds, out_path, **kwargs):
    """Reprojects a raster dataset to a new spatial reference system.

    Args:
        ds (gdal.Dataset or str): The input raster dataset or its path.
        out_path (str): The path to the output reprojected raster file.
        **kwargs: Additional arguments for gdal.WarpOptions().

    Returns:
        str: The path to the output reprojected raster file.
    """
    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    # input SpatialReference
    in_srs = kwargs.pop('srcSRS', WGS84)
    inSpatialRef = read_srs([ds, in_srs])

    # output SpatialReference
    out_srs = kwargs.pop('dstSRS', WGS84)
    outSpatialRef = read_srs(out_srs)

    resample_alg = kwargs.pop('resampleAlg', gdal.GRA_Average)
    warp_options = kwargs.pop('warpOptions', []) + ['UNIFIED_SRC_NODATA=PARTIAL']
    option = gdal.WarpOptions(
        creationOptions=CREATION,
        resampleAlg=resample_alg,
        srcSRS=inSpatialRef,
        dstSRS=outSpatialRef,
        multithread=True,
        warpOptions=warp_options,
        **kwargs,
    )
    gdal.Warp(out_file, ds, options=option)

    return out_file


def grib_to_tif(ds, out_path=None, **kwargs):
    """Converts a GRIB file to a GeoTIFF file.

    Args:
        ds (gdal.Dataset or str): The input GRIB dataset or its path.
        out_path (str, optional): The path to the output GeoTIFF file. If None, the output
                                  file is saved in the same directory as the input file.
                                  Defaults to None.
        **kwargs: Additional arguments for gdal.WarpOptions().

    Returns:
        str: The path to the output GeoTIFF file.
    """
    ds, ras = ds_name(ds)

    if os.path.splitext(os.path.basename(ras))[1] != '.grib':
        return

    if out_path:
        out_file = context_file(ras, out_path)
    else:
        out_file = os.path.join(
            os.path.dirname(ras), os.path.splitext(os.path.basename(ras))[0] + '.tif'
        )

    if os.path.exists(out_file):
        return out_file

    srcSRS = read_srs([ds, kwargs.pop('srcSRS', WGS84)])
    warp_options = kwargs.pop('warpOptions', []) + ['UNIFIED_SRC_NODATA=PARTIAL']
    option = gdal.WarpOptions(
        multithread=True,
        srcSRS=srcSRS,
        creationOptions=CREATION,
        warpOptions=warp_options,
        **kwargs,
    )
    gdal.Warp(out_file, ds, options=option)

    return out_file


def tif_copy_assign(out_file, ds_eg, array, srs=None, nodata=None):
    """Creates a new GeoTIFF file with the same georeferencing as an example dataset.

    Args:
        out_file (str): The path to the output GeoTIFF file.
        ds_eg (gdal.Dataset or str): The example raster dataset or its path.
        array (np.ndarray): The numpy array to write to the new file.
        srs (osr.SpatialReference, optional): The spatial reference system. If None, it is
                                              read from the example dataset. Defaults to None.
        nodata (float, optional): The nodata value. If None, it is read from the example
                                  dataset. Defaults to None.

    Returns:
        str: The path to the output GeoTIFF file.
    """
    out_file = str(out_file)
    if os.path.exists(out_file):
        return out_file
    ds_eg = ds_name(ds_eg)[0]

    if array.ndim == 2:
        array = array.reshape([1, *array.shape])
    if array.ndim != 3:
        raise (Exception('array dims must be 2 or 3'))

    # set nodata value
    if nodata is None:
        if ds_eg.GetRasterBand(1).GetNoDataValue() is not None:
            nodata = ds_eg.GetRasterBand(1).GetNoDataValue()
        else:
            raise (Exception('nodata must be passed'))
    if isinstance(array, np.ma.MaskedArray):
        array.set_fill_value(nodata)
        array = array.filled()

    ds = gdal.GetDriverByName('GTiff').Create(
        out_file,
        array.shape[2],
        array.shape[1],
        array.shape[0],
        TYPE_MAP[array.dtype.name],
        CREATION,
    )

    # fill with array
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    ds.WriteArray(array)

    # set geotransform
    trans = ds_eg.GetGeoTransform()
    ds.SetGeoTransform(tuple(trans))

    # set SpatialReference
    ds.SetSpatialRef(read_srs([srs, ds_eg]))

    return out_file
