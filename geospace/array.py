import numpy as np
from osgeo import gdal
from osgeo import gdal_array


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


def fill_nodata(arr, valid, nodata=np.nan):
    null = -9999 if np.isnan(nodata) else nodata
    arr[np.isnan(arr)] = null
    invalid = np.broadcast_to(~valid.astype(bool), arr.shape)
    mask = ((arr != null) | invalid).astype(np.uint8)
    dst_ds = gdal_array.OpenArray(arr)
    mask_ds = gdal_array.OpenArray(mask)
    for i in range(dst_ds.RasterCount):
        band = dst_ds.GetRasterBand(i + 1)
        maskBand = mask_ds.GetRasterBand(i + 1)
        gdal.FillNodata(
            targetBand=band,
            maskBand=maskBand,
            maxSearchDist=100,
            smoothingIterations=0,
            options={"NODATA": null},
        )
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
