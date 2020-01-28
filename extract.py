import os
import numpy as np
from osgeo import gdal, ogr, osr
from ._const import CREATION, CONFIG
from .transform import geo2imagexy, rep_file
from .block_write import block_write
from .create_tif import zeros_tif
from .projection import proj_shapefile


__all__ = ['extract']


# Now Read the large raster block by block, expecting to see no increase of
# memory as we loop through the blocks.
def map_burn(in_datas):
    in_datas[1].set_fill_value(0)
    burn_data = in_datas[1].filled()
    rect_data = np.ma.masked_where(np.logical_not(
        burn_data.astype(np.bool)), in_datas[0])
    return rect_data.filled()


def clip(ds, outLayer, no_data=None, rect_file=None,
         enlarge=10, save_cache=False, ext=None, new=True):

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
                            1], options=['ALL_TOUCHED=TRUE'])
        poly_ds = None

        option = gdal.WarpOptions(multithread=True, options=CONFIG,
                                  creationOptions=CREATION,
                                  xRes=rect.GetGeoTransform()[1],
                                  yRes=rect.GetGeoTransform()[5],
                                  resampleAlg=gdal.GRA_Average,
                                  outputType=gdal.GDT_Float32)
        burn_ds = gdal.Warp(burn_file, poly_file, options=option)
        burn_ds.GetRasterBand(1).SetNoDataValue(0)

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
        block_write(rect, [rect_band, burn_band], rect_band, map_burn)

    burn_ds = None
    poly_ds = None

    return rect, burn_data


def extract_stat(ras, shp, PROJ=None, no_data=None, **kwargs):
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
        # get the output layer's feature definition
        outLayerDefn = outLayer.GetLayerDefn()

        # add fields
        for i in range(0, inLayerDefn.GetFieldCount()):
            fieldDefn = inLayerDefn.GetFieldDefn(i)
            outLayer.CreateField(fieldDefn)

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
            stat[r, c - 1] = np.average(rect_data[select],
                                        weights=burn_data[select])

        outFeature = None
        outDataSet = None
        outLayer = None
        inFeature = None

    ds = None
    return stat


def extract(ras, shp, PROJ=None, no_data=None, stat=False, **kwargs):
    if stat:
        return extract_stat(ras, shp, PROJ=None, no_data=None, **kwargs)

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
    proj_shapefile(shp, out_shp, out_proj=outSpatialRef)
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
