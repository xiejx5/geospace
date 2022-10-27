from osgeo import gdal

# gdal config
gdal.PushErrorHandler('CPLQuietErrorHandler')
CREATION = ['BIGTIFF=YES', 'TILED=YES', 'NUM_THREADS=ALL_CPUS',
            'COMPRESS=ZSTD', 'PREDICTOR=1', 'ZSTD_LEVEL=1']

# mapping between gdal type and ogr field type
TYPE_MAP = {'uint8': gdal.GDT_Byte,
            'int8': gdal.GDT_Byte,
            'uint16': gdal.GDT_UInt16,
            'int16': gdal.GDT_Int16,
            'uint32': gdal.GDT_UInt32,
            'int32': gdal.GDT_Int32,
            'float32': gdal.GDT_Float32,
            'float64': gdal.GDT_Float64,
            'complex64': gdal.GDT_CFloat32,
            'complex128': gdal.GDT_CFloat64}
