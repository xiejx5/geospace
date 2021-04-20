from osgeo import gdal, osr
from ._const import CREATION


def zeros_tif(out_file, x_size, y_size, n_band,
              data_type, trans, srs, no_data=2):

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
    ds.SetProjection(inSpatialRef.ExportToWkt())

    band = None
    ds = None

    return out_file
