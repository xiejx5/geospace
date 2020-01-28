from osgeo import gdal, osr
from ._const import CREATION
from .transform import ds_name, context_file


def zeros_tif(out_path, x_size, y_size, n_band,
              data_type, trans, srs, no_data=2):
    ds, ras = ds_name(out_path)
    out_file = context_file(ras, out_path)

    ds = gdal.GetDriverByName('GTiff').Create(
        out_file, x_size, y_size, n_band, data_type, CREATION)

    # fill with 0 and set no data
    band = ds.GetRasterBand(1)
    band.Fill(0)
    band.SetNoDataValue(no_data)

    # set geotransform
    ds.SetGeoTransform(tuple(trans))

    # set SpatialReference
    if isinstance(srs, osr.SpatialReference):
        inSpatialRef = srs
    else:
        inSpatialRef = osr.SpatialReference()
        inSpatialRef.ImportFromProj4(srs)
    ds.SetProjection(srs.ExportToWkt())

    band = None
    ds = None

    return out_file
