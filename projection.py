import os
from osgeo import ogr, osr

__all__ = ['proj_shapefile']


def proj_shapefile(in_shp, out_shp, in_proj=None,
                   out_proj="+proj=longlat +datum=WGS84 +ellps=WGS84"):
    # Filename of input OGR file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if isinstance(in_shp, ogr.Layer):
        inLayer = in_shp
    else:
        source_ds = driver.Open(in_shp)
        inLayer = source_ds.GetLayer()

    # input SpatialReference
    inSpatialRef = inLayer.GetSpatialRef()
    if not inSpatialRef:
        if in_proj is not None:
            if isinstance(in_proj, osr.SpatialReference):
                inSpatialRef = in_proj
            else:
                inSpatialRef = osr.SpatialReference()
                inSpatialRef.ImportFromProj4(in_proj)
        else:
            raise(ValueError("in_proj must be set"))

    # output SpatialReference
    if isinstance(out_proj, osr.SpatialReference):
        outSpatialRef = out_proj
    elif os.path.isfile(out_proj):
        ds = driver.Open(out_proj)
        outSpatialRef = ds.GetLayer().GetSpatialRef()
        ds = None
    else:
        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromProj4(out_proj)

    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

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

    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
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
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()

    source_ds = None
    inLayer = None
    outLayer = None
    outDataSet = None

    return out_shp
