import os
import numpy as np
from osgeo import gdal, ogr
from functools import partial
from geospace._const import CREATION
from multiprocessing import Pool, cpu_count
from geospace.projection import read_srs
from geospace.spatial_calc import real_area
from geospace.raster import tif_copy_assign
from geospace.boundary import _enlarge_bound
from geospace.shape import shp_projection, shp_filter
from geospace.utils import zeros_tif, block_write, ds_name

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


def _clip(ds, outLayer, rect_file=None, enlarge=10, save_cache=False,
          reuse_cache=False, ext='', rasterize_option=['ALL_TOUCHED=TRUE']):
    # get no data
    no_data = ds.GetRasterBand(1).GetNoDataValue()

    # get shp extent in form of raster grids
    t = ds.GetGeoTransform()
    x_min, x_max, y_min, y_max = outLayer.GetExtent()
    bound = _enlarge_bound(ds, x_min, y_min, x_max, y_max)

    # clip with rectangle, generate rect_file
    if rect_file is None:
        rect_file = '/vsimem/_rect.tif'
    option = gdal.WarpOptions(multithread=True, outputBounds=bound,
                              dstSRS=outLayer.GetSpatialRef(),
                              creationOptions=CREATION, dstNodata=no_data,
                              xRes=t[1], yRes=t[5], srcNodata=no_data,
                              resampleAlg=gdal.GRA_NearestNeighbour)
    rect = gdal.Warp(rect_file, ds, options=option)
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


def extract(ras, shp, ras_srs="+proj=longlat +datum=WGS84 +ellps=WGS84",
            no_data=None, stat=False, **kwargs):
    ds, _ = ds_name(ras)

    # set projection
    SpatialRef = read_srs([ds, ras_srs])
    ds.SetProjection(SpatialRef.ExportToWkt())
    out_shp = '/vsimem/_outline.shp'
    shp_projection(shp, out_shp, out_srs=SpatialRef)
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
    rect, burn_data = _clip(ds, outLayer, **kwargs)
    if not stat:
        return rect.GetDescription()

    # initialize output for statistics
    stat = np.full(ds.RasterCount, np.nan)

    # iterate all bands
    rect_data = rect.ReadAsArray()
    mask = rect_data == no_data
    arr = np.ma.masked_array(rect_data, mask)

    try:
        if arr.ndim > 2:
            return np.ma.average(arr.reshape(arr.shape[0], -1),
                                 weights=burn_data.ravel(),
                                 axis=1).filled(np.nan)
        else:
            return np.ma.average(arr, weights=burn_data)
    except ZeroDivisionError:
        return


def basin_average_worker(shp, rasters, s, t, field, filter, **kwargs):
    filter_sql = f"{field} = '{filter}'"
    filter_shp = shp_filter(shp, filter_sql)
    one_out = np.full(t[-1], np.nan)
    for i, ras in enumerate(rasters):
        one_out[s[i]:t[i]] = extract(ras, filter_shp, stat=True,
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

    n_bands = [gdal.Open(ras).RasterCount for ras in rasters]
    t = np.cumsum(n_bands)
    s = np.roll(t, 1)
    s[0] = 0

    with Pool(min(cpu_count() * 3 // 4, len(filter))) as p:
        output = p.map(partial(basin_average_worker, shp, rasters,
                               s, t, field, **kwargs), filter)

    names = np.zeros(t[-1], dtype='object')
    for i, ras in enumerate(rasters):
        string = os.path.splitext(os.path.basename(ras))[0]
        names[s[i]:t[i]] = np.core.defchararray.add(string, np.char.mod('%d', np.arange(n_bands[i])))
        names[s[i]] = string

    return pd.DataFrame(output, columns=names, index=filter)
