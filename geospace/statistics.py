import os
import shutil
import numpy as np
from functools import partial
from osgeo import gdal, ogr, osr
from geospace._const import CREATION, CONFIG
from multiprocessing import Pool, cpu_count
from geospace.shape import project_shape, shp_filter
from geospace.utils import geo2imagexy, rep_file, zeros_tif, block_write

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


def clip(ds, outLayer, no_data=None, rect_file=None, enlarge=10,
         save_cache=False, ext=None, new=True, rasterize_option=['ALL_TOUCHED=TRUE']):

    # Open the data source and read in the extent
    t = ds.GetGeoTransform()
    x_min, x_max, y_min, y_max = outLayer.GetExtent()
    ulX, ulY = geo2imagexy(ds, x_min, y_min)
    lrX, lrY = geo2imagexy(ds, x_max, y_max)
    clip_range = [min(ulX, lrX), min(ulY, lrY),
                  abs(ulX - lrX) + 1, abs(ulY - lrY) + 1]
    ul_lon = t[0] + t[1] * clip_range[0] + t[2] * clip_range[1]
    ul_lat = t[3] + t[4] * clip_range[0] + t[5] * clip_range[1]
    lr_lon = t[0] + t[1] * (clip_range[0] + clip_range[2]) + \
        t[2] * (clip_range[1] + clip_range[3])
    lr_lat = t[3] + t[4] * (clip_range[0] + clip_range[2]) + \
        t[5] * (clip_range[1] + clip_range[3])
    bound = [min(ul_lon, lr_lon), min(ul_lat, lr_lat),
             max(ul_lon, lr_lon), max(ul_lat, lr_lat)]

    if save_cache:
        cache_dir = 'cache'
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)

        if ext is None:
            poly_file = rep_file(cache_dir, 'poly.tif')
            burn_file = rep_file(cache_dir, 'burn.tif')
        elif new:
            poly_file = rep_file(cache_dir, str(ext) + '_poly.tif')
            burn_file = rep_file(cache_dir, str(ext) + '_burn.tif')
        else:
            poly_file = os.path.join(cache_dir, str(ext) + '_poly.tif')
            burn_file = os.path.join(cache_dir, str(ext) + '_burn.tif')
    else:
        poly_file = '/vsimem/_poly.tif'
        burn_file = '/vsimem/_burn.tif'

    # set no data
    if ds.GetRasterBand(1).GetNoDataValue() is not None:
        no_data = ds.GetRasterBand(1).GetNoDataValue()
    if no_data is None:
        raise(ValueError("no_data must be initialed"))

    # create temp bool in_poly tif
    has_old = os.path.exists(burn_file) and not new

    # clip with rectangle
    if rect_file is None:
        rect_file = '/vsimem/_rect.tif'
    option = gdal.WarpOptions(multithread=True, options=CONFIG,
                              creationOptions=CREATION, outputBounds=bound,
                              dstSRS=outLayer.GetSpatialRef(), dstNodata=no_data,
                              xRes=t[1], yRes=t[5], srcNodata=no_data,
                              resampleAlg=gdal.GRA_NearestNeighbour)
    rect = gdal.Warp(rect_file, ds, options=option)

    if not has_old:
        # set geotransform
        trans = list(rect.GetGeoTransform())
        trans[1] = trans[1] / enlarge
        trans[5] = trans[5] / enlarge

        # set SpatialReference
        srs = outLayer.GetSpatialRef()

        zeros_tif(poly_file, int(clip_range[2] * enlarge),
                  int(clip_range[3] * enlarge), 1,
                  gdal.GDT_Byte, trans, srs, no_data=2)
        poly_ds = gdal.Open(poly_file, gdal.GA_Update)

        # Rasterize
        gdal.RasterizeLayer(poly_ds, [1], outLayer, burn_values=[
                            1], options=rasterize_option)
        poly_ds = None

        option = gdal.WarpOptions(multithread=True, options=CONFIG,
                                  creationOptions=CREATION, dstNodata=0,
                                  xRes=rect.GetGeoTransform()[1],
                                  yRes=rect.GetGeoTransform()[5],
                                  resampleAlg=gdal.GRA_Average,
                                  outputType=gdal.GDT_Float32)
        burn_ds = gdal.Warp(burn_file, poly_file, options=option)

        # return bool matrix in polygon
        burn_band = burn_ds.GetRasterBand(1)
        burn_data = burn_band.ReadAsArray()
    else:
        burn_ds = gdal.Open(burn_file)
        burn_band = burn_ds.GetRasterBand(1)
        burn_data = burn_band.ReadAsArray()

    # change rect
    for c in range(1, rect.RasterCount + 1):
        rect_band = rect.GetRasterBand(c)
        if rect_band.GetNoDataValue() is None:
            rect_band.SetNoDataValue(no_data)
        block_write(rect, [rect_band, burn_band], rect_band, _map_burn)

    burn_ds = None
    poly_ds = None

    return rect, burn_data


def _extract_stat(ras, shp, PROJ=None, no_data=None, **kwargs):
    if isinstance(ras, str):
        ds = gdal.Open(ras)
    else:
        ds = ras

    # Filename of input OGR file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    source_ds = driver.Open(shp)
    inLayer = source_ds.GetLayer()
    inLayerDefn = inLayer.GetLayerDefn()

    # input SpatialReference
    inSpatialRef = inLayer.GetSpatialRef()

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    if ds.GetProjectionRef():
        prj_srs = ds.GetProjectionRef()
        outSpatialRef.ImportFromWkt(prj_srs)
    elif PROJ:
        outSpatialRef.ImportFromProj4(PROJ)
        ds.SetProjection(outSpatialRef.ExportToWkt())
    else:
        raise(ValueError("PROJ must be initialed"))

    # create the CoordinateTransformation
    if int(gdal.__version__[0]) >= 3:
        outSpatialRef.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # initial output
    row = inLayer.GetFeatureCount()
    col = ds.RasterCount
    stat = np.full((row, col), np.nan)

    for r in range(row):
        inFeature = inLayer.GetFeature(r)
        # create the output layer
        outFile = '/vsimem/outline.shp'
        if os.path.exists(outFile):
            driver.DeleteDataSource(outFile)
        outDataSet = driver.CreateDataSource(outFile)
        outLayer = outDataSet.CreateLayer(outFile, outSpatialRef)

        # add fields
        for i in range(0, inLayerDefn.GetFieldCount()):
            fieldDefn = inLayerDefn.GetFieldDefn(i)
            outLayer.CreateField(fieldDefn)

        # get the output layer's feature definition
        outLayerDefn = outLayer.GetLayerDefn()

        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(
                i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)

        if ds.GetRasterBand(1).GetNoDataValue() is not None:
            no_data = ds.GetRasterBand(1).GetNoDataValue()
        if no_data is None:
            raise(ValueError("no_data must be initialed"))
        rect = None
        rect, burn_data = clip(ds, outLayer, no_data=no_data, **kwargs)

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

        outFeature = None
        outDataSet = None
        outLayer = None
        inFeature = None

    ds = None
    return stat


def extract(ras, shp, PROJ=None, no_data=None, stat=False, **kwargs):
    if stat:
        return _extract_stat(ras, shp, PROJ=PROJ, no_data=no_data, **kwargs)

    if isinstance(ras, str):
        ds = gdal.Open(ras)
    else:
        ds = ras

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    if ds.GetProjectionRef():
        prj_srs = ds.GetProjectionRef()
        outSpatialRef.ImportFromWkt(prj_srs)
    elif PROJ:
        prj_srs = PROJ
        outSpatialRef.ImportFromProj4(prj_srs)
        ds.SetProjection(outSpatialRef.ExportToWkt())
    else:
        raise(ValueError("PROJ must be initialed"))

    out_shp = '/vsimem/outline.shp'
    project_shape(shp, out_shp, out_srs=outSpatialRef)
    outDataSet = ogr.Open(out_shp)
    outLayer = outDataSet.GetLayer()

    # set no data
    if ds.GetRasterBand(1).GetNoDataValue() is not None:
        no_data = ds.GetRasterBand(1).GetNoDataValue()
    if no_data is None:
        raise(ValueError("no_data must be initialed"))

    # clip with shapefile
    rect, _ = clip(ds, outLayer, no_data=no_data, **kwargs)

    return rect.GetDescription()


def shp_weighted_mean(in_shp, clip_shp, field, out_shp=None, save_cache=False):
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # get layer of in_shp
    if isinstance(in_shp, str):
        ds = driver.Open(in_shp)
    else:
        ds = in_shp

    in_layer = ds.GetLayer()
    srs = ds.GetLayer().GetSpatialRef()

    # project clip_shp
    if save_cache:
        proj_shp = rep_file('cache', os.path.splitext(
            os.path.basename(clip_shp))[0] + '_proj.shp')
    else:
        proj_shp = '/vsimem/_proj.shp'
    project_shape(clip_shp, proj_shp, out_srs=srs)
    clip_ds = driver.Open(proj_shp)
    clip_layer = clip_ds.GetLayer()

    # export out_shp
    if out_shp is None:
        if save_cache:
            out_shp = rep_file('cache', os.path.splitext(
                os.path.basename(clip_shp))[0] + '_out.shp')
        else:
            out_shp = '/vsimem/_out.shp'

    out_ds = driver.CreateDataSource(out_shp)
    out_layer = out_ds.CreateLayer(out_shp, srs=srs)

    in_layer.Clip(clip_layer, out_layer)

    area = []
    logK = []
    # newField = ogr.FieldDefn('Area', ogr.OFTReal)
    # out_layer.CreateField(newField)
    c = out_layer.GetFeatureCount()
    for i in range(c):
        f = out_layer.GetFeature(i)
        area.append(f.GetGeometryRef().GetArea())
        logK.append(f.GetField(field))
        # f.SetField('Area', f.GetGeometryRef().GetArea())
        # out_layer.SetFeature(f)
    area = np.array(area)
    logK = np.array(logK)
    mean_logK = np.average(logK, weights=area)
    # mean_logK = np.log10(np.average(np.power(10, logK / 100), weights=area))
    out_layer = None
    out_ds = None
    return mean_logK


def basin_average_worker(shp, rasters, basin_id, **kwargs):
    filter_sql = f"STAID = '{basin_id}'"
    filter_shp = shp_filter(shp, filter_sql)
    one_out = np.full(len(rasters), np.nan)
    for i, ras in enumerate(rasters):
        one_out[i] = np.squeeze(extract(ras, filter_shp, enlarge=10,
                                        ext=basin_id, stat=True, **kwargs))
    return one_out


def basin_average(shp, rasters, basins_id, **kwargs):
    if isinstance(rasters, str):
        rasters = [rasters]
    with Pool(cpu_count() * 3 // 4) as p:
        output = p.map(partial(basin_average_worker, shp, rasters, **kwargs), basins_id)
    if kwargs.pop('save_cache', False):
        if os.path.exists('cache'):
            shutil.rmtree('cache')

    names = [os.path.splitext(os.path.basename(i))[0] for i in rasters]
    return pd.DataFrame(output, columns=names)
