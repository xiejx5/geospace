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
    srcSRS = kwargs.pop("srcSRS", "EPSG:4326")
    dstSRS = kwargs.pop("dstSRS", "EPSG:4326")
    resampleAlg = kwargs.pop("resampleAlg", gdal.GRA_Average)

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
    if isinstance(ds, gdal.Dataset):
        trans, crs = ds.GetGeoTransform(), ds.GetProjection()
        arr = ds.ReadAsArray()
    else:
        trans, crs = shape_to_trans(*ds.shape[-2:]), "EPSG:4326"
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

    kwargs.setdefault("maxSearchDist", 50)
    kwargs.setdefault("smoothingIterations", 0)
    kwargs.setdefault("options", {"NODATA": null, "TEMP_FILE_DRIVER": "MEM"})
    kwargs["options"]["INTERPOLATION"] = "NEAREST" if nearest else "INV_DIST"

    for i in range(dst_ds.RasterCount):
        band = dst_ds.GetRasterBand(i + 1)
        band.SetNoDataValue(null)
        maskBand = mask_ds.GetRasterBand(i + 1)
        gdal.FillNodata(targetBand=band, maskBand=maskBand, **kwargs)

    if (arr[~(mask | invalid)] == null).any():
        if not fast:
            # Some nodata not filled, consider larger maxSearchDist
            kwargs["maxSearchDist"] *= 2
        arr = fill_nodata(ds, valid, nodata, nearest, fast=False, **kwargs)
    arr[invalid] = nodata
    return arr


if __name__ == "__main__":
    # Example usage
    arr = reproject(
        np.zeros((3600, 7200), dtype=np.float32),
        (1800, 3600),
        (-180, 0.05, 0, 90, 0, -0.05),
        (-180, 0.1, 0, 90, 0, -0.1),
    )
