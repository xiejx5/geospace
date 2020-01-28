import os
import locale
from osgeo import gdal, ogr, osr


__all__ = ['split_shp']


def split_shp(outDir, shp_paths):
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
    gdal.SetConfigOption("SHAPE_ENCODING", locale.getpreferredencoding())
    for shp in shp_paths:
        # Filename of input OGR file
        driver = ogr.GetDriverByName("ESRI Shapefile")
        source_ds = driver.Open(shp)
        inLayer = source_ds.GetLayer()

        # input SpatialReference
        inSpatialRef = inLayer.GetSpatialRef()

        # output SpatialReference
        outSpatialRef = osr.SpatialReference()
        prj_srs = "+proj=longlat +datum=WGS84 +ellps=WGS84"
        outSpatialRef.ImportFromProj4(prj_srs)

        # create the CoordinateTransformation
        coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

        # loop through the input features
        inLayerDefn = inLayer.GetLayerDefn()
        inFeature = inLayer.GetNextFeature()
        while inFeature:
            # create the output layer
            outFile = outDir + '/' + inFeature.GetField(0) + '.shp'
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
            # dereference the features and get the next input feature
            outFeature = None
            outDataSet = None
            outLayer = None
            inFeature = inLayer.GetNextFeature()

        outLayer = None
        inLayer = None
        source_ds = None
        outDataSet = None
