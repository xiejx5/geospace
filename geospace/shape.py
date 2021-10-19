import os
import locale
from osgeo import gdal, ogr, osr


def shp_buffer(in_shp, out_shp, buffdist, in_srs=None):
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
        if in_srs is not None:
            if isinstance(in_srs, osr.SpatialReference):
                inSpatialRef = in_srs
            else:
                inSpatialRef = osr.SpatialReference()
                inSpatialRef.ImportFromProj4(in_srs)
        else:
            raise(ValueError("in_srs must be set"))

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


def proj_shapefile(in_shp, out_shp, in_srs=None,
                   out_srs="+proj=longlat +datum=WGS84 +ellps=WGS84"):
    gdal.SetConfigOption("SHAPE_ENCODING", 'utf-8')

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
        if in_srs is not None:
            if isinstance(in_srs, osr.SpatialReference):
                inSpatialRef = in_srs
            else:
                inSpatialRef = osr.SpatialReference()
                inSpatialRef.ImportFromProj4(in_srs)
        else:
            raise(ValueError("in_srs must be set"))

    # output SpatialReference
    if isinstance(out_srs, osr.SpatialReference):
        outSpatialRef = out_srs
    elif os.path.isfile(out_srs):
        ds = driver.Open(out_srs)
        outSpatialRef = ds.GetLayer().GetSpatialRef()
        ds = None
    else:
        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromProj4(out_srs)

    # create the CoordinateTransformation
    if int(gdal.__version__[0]) >= 3:
        outSpatialRef.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
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


def split_shp(outDir, shp_paths, name_index=1):
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
    gdal.SetConfigOption("SHAPE_ENCODING", locale.getpreferredencoding())
    if not isinstance(shp_paths, list):
        shp_paths = [shp_paths]
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
        if int(gdal.__version__[0]) >= 3:
            outSpatialRef.SetAxisMappingStrategy(
                osr.OAMS_TRADITIONAL_GIS_ORDER)
        coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

        # loop through the input features
        inLayerDefn = inLayer.GetLayerDefn()
        inFeature = inLayer.GetNextFeature()
        while inFeature:
            # create the output layer
            outFile = outDir + '/' + inFeature.GetField(name_index) + '.shp'
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


def shp_filter(shps, filter_sql, filter_shp=None):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds_shp = driver.Open(shps, 0)
    layer = ds_shp.GetLayer()
    layer.SetAttributeFilter(filter_sql)
    if filter_shp is None:
        filter_shp = '/vsimem/filter.shp'
    filter = driver.CreateDataSource(filter_shp)
    filter.CopyLayer(layer, 'filter')
    filter = None
    return filter_shp
