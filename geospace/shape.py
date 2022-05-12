import os
import re
from osgeo import gdal, ogr
from geospace.projection import read_srs, coord_trans


def shp_buffer(in_shp, out_shp, buffdist, in_srs=None):
    return shp_geom_map(in_shp, out_shp, in_srs=in_srs,
                        func=lambda geom: geom.Buffer(float(buffdist)))


def shp_projection(in_shp, out_shp, in_srs=None,
                   out_srs="+proj=longlat +datum=WGS84 +ellps=WGS84"):
    gdal.SetConfigOption("SHAPE_ENCODING", 'utf-8')

    # Filename of input OGR file
    if isinstance(in_shp, ogr.Layer):
        inLayer = in_shp
    else:
        source_ds = ogr.Open(in_shp)
        inLayer = source_ds.GetLayer()

    # input SpatialReference
    inSpatialRef = read_srs([inLayer, in_srs])

    # output SpatialReference
    outSpatialRef = read_srs(out_srs)

    # create the CoordinateTransformation
    coordTrans = coord_trans(inSpatialRef, outSpatialRef)

    return shp_geom_map(in_shp, out_shp,
                        out_srs=outSpatialRef,
                        func=lambda geom: geom.Transform(coordTrans))


def shp_filter(shps, filter_sql, filter_shp=None):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds_shp = ogr.Open(shps, 0)
    layer = ds_shp.GetLayer()
    if filter_shp is None:
        filter_shp = '/vsimem/filter.shp'

    res = re.findall(r"None = '([0-9]+)'", filter_sql)
    if len(res) > 0:
        return shp_geom_map(layer, filter_shp, idxs=int(res[0]))

    layer.SetAttributeFilter(filter_sql)
    filter = driver.CreateDataSource(filter_shp)
    filter.CopyLayer(layer, 'filter')
    return filter_shp


def shp_geom_map(in_shp, out_shp, in_srs=None,
                 out_srs=None, idxs=None, func=None):
    gdal.SetConfigOption("SHAPE_ENCODING", 'utf-8')

    # Filename of input OGR file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if isinstance(in_shp, ogr.Layer):
        inLayer = in_shp
    else:
        source_ds = ogr.Open(in_shp)
        inLayer = source_ds.GetLayer()

    # output SpatialReference
    outSpatialRef = read_srs([out_srs, inLayer, in_srs])

    # create the output layer
    if 'vsimem' not in os.path.dirname(out_shp):
        if (not os.path.exists(os.path.dirname(out_shp)) and
                os.path.dirname(out_shp) != ''):
            os.makedirs(os.path.dirname(out_shp))
    if os.path.exists(out_shp):
        driver.DeleteDataSource(out_shp)
    outDataSet = driver.CreateDataSource(out_shp)
    outLayer = outDataSet.CreateLayer(out_shp, outSpatialRef)

    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()

    # set all feature idxs if not assigned
    if idxs is None:
        idxs = range(inLayer.GetFeatureCount())
    elif isinstance(idxs, int):
        idxs = [idxs]

    # loop through the input features
    for i in idxs:
        inFeature = inLayer.GetFeature(i)
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom if func is None else func(geom))
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(
                i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None

    return out_shp
