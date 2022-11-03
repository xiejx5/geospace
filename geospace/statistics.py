import os
import numpy as np
from osgeo import gdal, ogr
from functools import partial
from geospace._const import WGS84, CREATION
from multiprocessing import Pool, cpu_count
from geospace.projection import read_srs
from geospace.spatial_calc import real_area
from geospace.raster import tif_copy_assign
from geospace.boundary import _enlarge_bound
from geospace.shape import shp_projection, shp_filter
from geospace.utils import (zeros_tif, block_write,
                            context_file, ds_name, rep_name)

try:
    import pandas as pd
except Exception:
    pass


# Now Read the large raster block by block, expecting to see no increase of
# memory as we loop through the blocks.
def _map_burn(arrs, burn_data):
    mask = np.broadcast_to(np.logical_not(burn_data.astype(bool)), arrs[0].shape)
    rect_data = np.ma.masked_where(mask, arrs[0])
    return rect_data.filled()


def _clip(ds, outLayer, out_file, enlarge=10, ext='', save_cache=False,
          reuse_cache=False, rasterize_option=['ALL_TOUCHED=TRUE']):
    # get no data
    no_data = ds.GetRasterBand(1).GetNoDataValue()

    # get shp extent in form of raster grids
    t = ds.GetGeoTransform()
    x_min, x_max, y_min, y_max = outLayer.GetExtent()
    bound = _enlarge_bound(ds, x_min, y_min, x_max, y_max)

    # clip with rectangle out_file, use Warp instead of Translate
    # to project longitude from (0, 360) to (-180, 180)
    option = gdal.WarpOptions(multithread=True, outputBounds=bound,
                              srcSRS=outLayer.GetSpatialRef(),
                              dstSRS=outLayer.GetSpatialRef(),
                              creationOptions=CREATION, dstNodata=no_data,
                              xRes=t[1], yRes=t[5], srcNodata=no_data,
                              resampleAlg=gdal.GRA_NearestNeighbour)
    rect = gdal.Warp(out_file, ds, options=option)
    x_size, y_size = rect.RasterXSize, rect.RasterYSize

    # paths of cached poly.tif and burn.tif
    if save_cache:
        reuse_cache = True
        cache_dir = 'cache'
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
        burn_file = os.path.join(cache_dir, str(ext) + '_burn.tif')
    elif reuse_cache:
        burn_file = os.path.join('/vsimem/', str(ext) + '_burn.tif')
    else:
        burn_file = '/vsimem/_burn_renew.tif'

    # create boolean poly.tif and burn.tif
    burn_ds = gdal.Open(burn_file)
    if reuse_cache and (burn_ds is not None):
        burn_data = burn_ds.ReadAsArray()
    else:
        # set geotransform
        trans = list(rect.GetGeoTransform())
        trans[1] = trans[1] / enlarge
        trans[5] = trans[5] / enlarge

        # Rasterize
        poly_file = '/vsimem/_poly.tif'
        zeros_tif(poly_file, x_size * enlarge, y_size * enlarge, 1,
                  gdal.GDT_Byte, trans, outLayer.GetSpatialRef(), no_data=2)
        poly_ds = gdal.Open(poly_file, gdal.GA_Update)
        gdal.RasterizeLayer(poly_ds, [1], outLayer, burn_values=[1],
                            options=rasterize_option)
        poly_data = poly_ds.ReadAsArray()

        # https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-in-python-using-numpy-d1bf0dd7b6f7
        if enlarge == 1:
            burn_data = poly_data
        else:
            burn_data = poly_data.reshape(y_size, enlarge,
                                          x_size, enlarge).swapaxes(1, 2)
            burn_data = burn_data.mean(axis=(-2, -1), dtype=np.float32)
            row_area = real_area(rect, np.arange(y_size), return_row_area=True)
            burn_data = burn_data * row_area.reshape([y_size, -1])

        # calculate shapefile intersection area in each grid
        if reuse_cache:
            tif_copy_assign(burn_file, rect, burn_data, no_data=0)

    # change value in the rectangle
    block_write(rect, _map_burn, burn_data)

    return rect, burn_data


def extract(ras, shp, out_path=None,
            ras_srs=WGS84, no_data=None, **kwargs):
    ds, ras = ds_name(ras)
    if out_path is None:
        out_file = '/vsimem/_rect.tif'
    else:
        out_file = context_file(ds.GetDescription(), out_path)
        if os.path.exists(out_file):
            return out_file
        kwargs['enlarge'] = 1

    # set projection
    out_shp = '/vsimem/_outline.shp'
    shp_projection(shp, out_shp, out_srs=read_srs([ds, ras_srs]))
    outDataSet = ogr.Open(out_shp)
    outLayer = outDataSet.GetLayer()

    # set no data
    if ds.GetRasterBand(1).GetNoDataValue() is not None:
        no_data = ds.GetRasterBand(1).GetNoDataValue()
    elif no_data is not None:
        ds.GetRasterBand(1).SetNoDataValue(no_data)
    else:
        raise (ValueError("no_data must be initialized"))

    # clip with the whole shapefile
    rect, burn_data = _clip(ds, outLayer, out_file, **kwargs)
    if out_path is not None:
        return rect.GetDescription()

    # iterate all bands
    rect_data = rect.ReadAsArray()
    mask = rect_data == no_data
    arr = np.ma.masked_array(rect_data, mask)
    arr = arr[np.newaxis, :, :] if arr.ndim != 3 else arr

    return np.ma.average(arr.reshape(arr.shape[0], -1),
                         weights=burn_data.ravel(),
                         axis=1).filled(np.nan)


def basin_average_worker(shp, rasters, s, t, field, filter, **kwargs):
    filter_sql = f"{field} = '{filter}'"
    filter_shp = shp_filter(shp, filter_sql)
    one_out = np.full(t[-1], np.nan)
    for i, ras in enumerate(rasters):
        one_out[s[i]:t[i]] = extract(ras, filter_shp,
                                     ext=filter, **kwargs)
    return one_out


def basin_average(shp, rasters, field='STAID', filter=None, **kwargs):
    if isinstance(rasters, str):
        rasters = [rasters]
    if filter is None:
        ds = ogr.Open(shp)
        layer = ds.GetLayer()
        filter = range(layer.GetFeatureCount())
        field = None
    if isinstance(filter, str) or isinstance(filter, int):
        filter = [filter]

    names, s, t = rep_name(rasters)

    with Pool(min(cpu_count() * 3 // 4, len(filter))) as p:
        output = p.map(partial(basin_average_worker, shp, rasters,
                               s, t, field, **kwargs), filter)

    return pd.DataFrame(output, columns=names, index=filter)
