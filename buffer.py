import os
from osgeo import ogr, osr

__all__ = ['shp_buffer']


def shp_buffer(in_shp, out_shp, buffdist, in_proj=None):
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

    # create the output layer
    if 'vsimem' not in os.path.dirname(out_shp):
        if (not os.path.exists(os.path.dirname(out_shp)) and
                os.path.dirname(out_shp) != ''):
            os.makedirs(os.path.dirname(out_shp))
    if os.path.exists(out_shp):
        driver.DeleteDataSource(out_shp)
    outDataSet = driver.CreateDataSource(out_shp)
    outLayer = outDataSet.CreateLayer(out_shp, inSpatialRef)

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
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom.Buffer(float(buffdist)))
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(
                i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()

    inLayer = None
    outLayer = None
    outDataSet = None
    source_ds = None

    return out_shp
