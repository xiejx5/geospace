import os
from osgeo import gdal, ogr
from ._const import CREATION, CONFIG


__all__ = ['shp_to_raster']


def shp_to_raster(shp, attr, out_path, ds_eg, tem_path, **kwargs):
    # create out put name
    out_file = os.path.join(out_path, os.path.splitext(
        os.path.basename(shp))[0] + '.tif')
    if os.path.exists(out_file):
        return
    tem_file = os.path.join(tem_path, os.path.splitext(
        os.path.basename(shp))[0] + '.tif')

    # extent warp options
    ds_ex = gdal.Translate('/vsimem/_extent.tif', ds_eg, bandList=[1])
    t = ds_eg.GetGeoTransform()
    temp_option = gdal.WarpOptions(multithread=True, options=CONFIG,
                                   creationOptions=CREATION, **kwargs,
                                   xRes=t[1] / 10, yRes=t[5] / 10,
                                   outputType=gdal.GDT_Float64)

    ds_tem = gdal.Warp(tem_file, ds_ex, options=temp_option)
    band = ds_tem.GetRasterBand(1)
    option = gdal.WarpOptions(multithread=True, options=CONFIG,
                              creationOptions=CREATION,  **kwargs,
                              xRes=t[1], yRes=t[5], resampleAlg=gdal.GRA_Average,
                              outputType=gdal.GDT_Float64)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp_factor = driver.Open(shp)
    layer = shp_factor.GetLayer()

    # create and use RasterizeLayer
    band.Fill(band.GetNoDataValue())
    gdal.RasterizeLayer(ds_tem, [1], layer,
                        options=["ATTRIBUTE=%s" % attr, 'ALL_TOUCHED=TRUE'])
    ds_out = gdal.Warp(out_file, ds_tem, options=option)

    # deal with units
    if os.path.splitext(os.path.basename(out_file))[0] == 'permeability':
        band_out = ds_out.GetRasterBand(1)
        no_data = band_out.GetNoDataValue()
        values = ds_out.ReadAsArray()
        values[values != no_data] = values[values != no_data] / 100
        band_out.WriteArray(values)

        band = ds_tem.GetRasterBand(1)
        no_data = band.GetNoDataValue()
        values = ds_tem.ReadAsArray()
        values[values != no_data] = values[values != no_data] / 100
        band.WriteArray(values)

    band_out = None
    ds_out = None
    band = None
    ds_tem = None
    shp_factor = None
    layer = None
