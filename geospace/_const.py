import os
from osgeo import gdal


# gdal config
gdal.DontUseExceptions()
gdal.SetConfigOption("SHAPE_ENCODING", 'utf-8')
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
gdal.SetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS")

# cpu used
N_CPU = max(int(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))
                * int(os.environ.get('SLURM_NTASKS', 1)) * 3 // 4), 1)

# default spatial reference system
WGS84 = "EPSG:4326"

# creation options
CREATION = ['BIGTIFF=YES', 'TILED=YES', 'NUM_THREADS=ALL_CPUS',
            'COMPRESS=ZSTD', 'PREDICTOR=1', 'ZSTD_LEVEL=1']

# mapping between gdal type and ogr field type
TYPE_MAP = {'bool': gdal.GDT_Byte,
            'uint8': gdal.GDT_Byte,
            'int8': gdal.GDT_Byte,
            'uint16': gdal.GDT_UInt16,
            'int16': gdal.GDT_Int16,
            'uint32': gdal.GDT_UInt32,
            'int32': gdal.GDT_Int32,
            'float32': gdal.GDT_Float32,
            'float64': gdal.GDT_Float64,
            'complex64': gdal.GDT_CFloat32,
            'complex128': gdal.GDT_CFloat64}
