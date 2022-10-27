import os
import shutil
import numpy as np
from osgeo import gdal, ogr
from functools import partial
from geospace._const import CREATION
from multiprocessing import Pool, cpu_count
from geospace.projection import read_srs
from geospace.spatial_calc import real_area
from geospace.boundary import _enlarge_bound
from geospace.shape import shp_projection, shp_filter, shp_geom_map
from geospace.utils import rep_file, zeros_tif, block_write, ds_name

try:
    import pandas as pd
except Exception:
    pass


# Now Read the large raster block by block, expecting to see no increase of
# memory as we loop through the blocks.
def _map_burn(in_datas):
    in_datas[1].set_fill_value(0)
    burn_data = in_datas[1].filled()
    rect_data = np.ma.masked_where(np.logical_not(
        burn_data.astype(np.bool)), in_datas[0])
    return rect_data.filled()


def _clip(ds, outLayer, rect_file=None, enlarge=10, save_cache=False,
          reuse_cache=False, ext='', rasterize_option=['ALL_TOUCHED=TRUE']):
    # get no data
    no_data = ds.GetRasterBand(1).GetNoDataValue()

    # get shp extent in form of raster grids
    t = ds.GetGeoTransform()
    x_min, x_max, y_min, y_max = outLayer.GetExtent()
    bound, clip_range = _enlarge_bound(ds, x_min, y_min, x_max, y_max)

    # clip with rectangle, generate rect_file
    if rect_file is None:
        rect_file = '/vsimem/_rect.tif'
    option = gdal.WarpOptions(multithread=True, outputBounds=bound,
                              dstSRS=outLayer.GetSpatialRef(),
                              creationOptions=CREATION, dstNodata=no_data,
                              xRes=t[1], yRes=t[5], srcNodata=no_data,
                              resampleAlg=gdal.GRA_NearestNeighbour)
    rect = gdal.Warp(rect_file, ds, options=option)

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
    if (not reuse_cache) or (burn_ds is None):
        # set geotransform
        trans = list(rect.GetGeoTransform())
        trans[1] = trans[1] / enlarge
        trans[5] = trans[5] / enlarge
        poly_file = '/vsimem/_poly.tif'
        zeros_tif(poly_file, int(clip_range[2] * enlarge),
                  int(clip_range[3] * enlarge), 1,
                  gdal.GDT_Byte, trans, outLayer.GetSpatialRef(), no_data=2)
        poly_ds = gdal.Open(poly_file, gdal.GA_Update)

        # Rasterize
        gdal.RasterizeLayer(poly_ds, [1], outLayer, burn_values=[
                            1], options=rasterize_option)
        poly_ds = None
        option = gdal.WarpOptions(multithread=True,
                                  creationOptions=CREATION, dstNodata=0,
                                  xRes=rect.GetGeoTransform()[1],
                                  yRes=rect.GetGeoTransform()[5],
                                  resampleAlg=gdal.GRA_Average,
                                  outputType=gdal.GDT_Float32)
        burn_ds = gdal.Warp(burn_file, poly_file, options=option)

    # return bool matrix in the shapefile polygon
    burn_band = burn_ds.GetRasterBand(1)
    burn_data = burn_band.ReadAsArray()

    # calculate real area for geodetic grids
    if 'PROJCS' not in burn_ds.GetProjection():
        row_area = real_area(burn_ds, np.arange(burn_data.shape[0]), return_row_area=True)
        burn_area = burn_data * row_area.reshape([burn_data.shape[0], -1])
        total_area = np.sum(burn_area)
        if total_area > 0:
            burn_data = burn_area / total_area

    # change value in the rectangle
    for c in range(1, rect.RasterCount + 1):
        rect_band = rect.GetRasterBand(c)
        if rect_band.GetNoDataValue() is None:
            rect_band.SetNoDataValue(no_data)
        block_write(rect, [rect_band, burn_band], rect_band, _map_burn)

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
    if not stat:
        return _clip(ds, outLayer, **kwargs)[0].GetDescription()

    # initialize output for statistics
    row, col = outLayer.GetFeatureCount(), ds.RasterCount
    stat = np.full((row, col), np.nan)

    for r in range(row):
        sin_shp = '/vsimem/_single.shp'
        sinDataSet = ogr.Open(shp_geom_map(outLayer, sin_shp, idxs=r))
        sinLayer = sinDataSet.GetLayer()
        rect, burn_data = _clip(ds, sinLayer, **kwargs)

        # iterate all bands
        for c in range(1, col + 1):
            rect_band = rect.GetRasterBand(c)
            rect_data = rect_band.ReadAsArray()
            no_data = rect_band.GetNoDataValue()
            select = (rect_data != no_data)
            try:
                stat[r, c - 1] = np.average(rect_data[select],
                                            weights=burn_data[select])
            except ZeroDivisionError:
                pass

    return stat


def basin_average_worker(shp, rasters, s, t, field, filter, **kwargs):
    filter_sql = f"{field} = '{filter}'"
    filter_shp = shp_filter(shp, filter_sql)
    one_out = np.full(t[-1], np.nan)
    for i, ras in enumerate(rasters):
        one_out[s[i]:t[i]] = np.squeeze(extract(ras, filter_shp, stat=True,
                                        ext=filter, **kwargs))
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
