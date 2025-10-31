import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from geospace.mask import shape_to_trans


def reproject(
    src_arr,
    dst_shape,
    src_transform,
    dst_transform,
    nodata=np.nan,
    **kwargs,
):
    """Reprojects a numpy array from a source to a destination projection.

    Args:
        src_arr (np.ndarray): The source numpy array to reproject.
        dst_shape (tuple): The shape of the destination array (rows, cols).
        src_transform (tuple): The geotransform of the source array.
        dst_transform (tuple): The geotransform of the destination array.
        nodata (float, optional): The nodata value. Defaults to np.nan.
        srcSRS (str, optional): The source spatial reference system. Defaults to 'EPSG:4326'.
        dstSRS (str, optional): The destination spatial reference system. Defaults to 'EPSG:4326'.
        resampleAlg (int, optional): The resampling algorithm. Defaults to gdal.GRA_Average.

    Returns:
        np.ndarray: The reprojected numpy array.
    """
    srcSRS = kwargs.pop('srcSRS', 'EPSG:4326')
    dstSRS = kwargs.pop('dstSRS', 'EPSG:4326')
    resampleAlg = kwargs.pop('resampleAlg', gdal.GRA_Average)

    option = gdal.WarpOptions(
        srcNodata=nodata,
        dstNodata=nodata,
        srcSRS=srcSRS,
        dstSRS=dstSRS,
        multithread=True,
        resampleAlg=resampleAlg,
        **kwargs,
    )

    src_ds = gdal_array.OpenArray(src_arr)
    src_ds.SetGeoTransform(src_transform)

    dst_arr = np.full(dst_shape, nodata, dtype=src_arr.dtype)
    dst_ds = gdal_array.OpenArray(dst_arr)
    dst_ds.SetGeoTransform(dst_transform)

    gdal.Warp(dst_ds, src_ds, options=option)
    return dst_arr


def fill_nodata(ds, valid, nodata=np.nan, nearest=False, fast=True, **kwargs):
    """Fills nodata values in a numpy array or GDAL dataset.

    This function uses GDAL's FillNodata algorithm to fill areas with no data.

    Args:
        ds (gdal.Dataset or np.ndarray): The input dataset or numpy array.
        valid (np.ndarray): A boolean array indicating valid data pixels.
        nodata (float, optional): The nodata value. Defaults to np.nan.
        nearest (bool, optional): If True, use nearest neighbor interpolation.
                                  If False, use inverse distance weighting. Defaults to False.
        fast (bool, optional): If True, performs a faster fill but might leave some holes.
                               If False, performs a more thorough fill. Defaults to True.
        maxSearchDist (int, optional): The maximum search distance. Defaults to 50.
        smoothingIterations (int, optional): The number of smoothing iterations. Defaults to 0.

    Returns:
        np.ndarray: The array with nodata values filled.
    """
    if isinstance(ds, gdal.Dataset):
        trans, crs = ds.GetGeoTransform(), ds.GetProjection()
        arr = ds.ReadAsArray()
    else:
        trans, crs = shape_to_trans(*ds.shape[-2:]), 'EPSG:4326'
        arr = ds.copy()

    null = -9999 if np.isnan(nodata) else nodata
    arr[np.isnan(arr)] = null
    dst_ds = gdal_array.OpenArray(arr)
    dst_ds.SetGeoTransform(trans)
    dst_ds.SetProjection(crs)

    invalid = np.broadcast_to(~valid.astype(bool), arr.shape)
    mask = ((arr != null) | invalid) if fast else (arr != null)
    mask_ds = gdal_array.OpenArray(mask.astype(np.uint8))
    mask_ds.SetGeoTransform(trans)
    mask_ds.SetProjection(crs)

    kwargs.setdefault('maxSearchDist', 50)
    kwargs.setdefault('smoothingIterations', 0)
    kwargs.setdefault('options', {'NODATA': null, 'TEMP_FILE_DRIVER': 'MEM'})
    kwargs['options']['INTERPOLATION'] = 'NEAREST' if nearest else 'INV_DIST'

    for i in range(dst_ds.RasterCount):
        band = dst_ds.GetRasterBand(i + 1)
        band.SetNoDataValue(null)
        maskBand = mask_ds.GetRasterBand(i + 1)
        gdal.FillNodata(targetBand=band, maskBand=maskBand, **kwargs)

    if (arr[~(mask | invalid)] == null).any():
        if not fast:
            # Some nodata not filled, consider larger maxSearchDist
            kwargs['maxSearchDist'] *= 2
        arr = fill_nodata(ds, valid, nodata, nearest, fast=False, **kwargs)
    arr[invalid] = nodata
    return arr


if __name__ == '__main__':
    # Example usage
    arr = reproject(
        np.zeros((3600, 7200), dtype=np.float32),
        (1800, 3600),
        (-180, 0.05, 0, 90, 0, -0.05),
        (-180, 0.1, 0, 90, 0, -0.1),
    )
