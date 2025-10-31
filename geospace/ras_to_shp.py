import os
from osgeo import gdal, ogr
from geospace.projection import read_srs


def polygonize(ds, shp_path):
    """Converts a raster dataset to a polygon shapefile.

    Args:
        ds (gdal.Dataset): The input raster dataset.
        shp_path (str): The path to the output shapefile.

    Returns:
        str: The path to the output shapefile.
    """
    if os.path.exists(shp_path):
        return
    # mapping between gdal type and ogr field type
    type_mapping = {
        gdal.GDT_Byte: ogr.OFTInteger,
        gdal.GDT_UInt16: ogr.OFTInteger,
        gdal.GDT_Int16: ogr.OFTInteger,
        gdal.GDT_UInt32: ogr.OFTInteger,
        gdal.GDT_Int32: ogr.OFTInteger,
        gdal.GDT_Float32: ogr.OFTReal,
        gdal.GDT_Float64: ogr.OFTReal,
        gdal.GDT_CInt16: ogr.OFTInteger,
        gdal.GDT_CInt32: ogr.OFTInteger,
        gdal.GDT_CFloat32: ogr.OFTReal,
        gdal.GDT_CFloat64: ogr.OFTReal,
    }

    srcband = ds.GetRasterBand(1)
    dst_layername = 'Shape'
    drv = ogr.GetDriverByName('ESRI Shapefile')
    dst_ds = drv.CreateDataSource(shp_path)

    dst_layer = dst_ds.CreateLayer(dst_layername, srs=read_srs(ds))
    raster_field = ogr.FieldDefn('id', type_mapping[srcband.DataType])
    dst_layer.CreateField(raster_field)
    gdal.Polygonize(srcband, srcband, dst_layer, 0, [], callback=None)
    return shp_path
