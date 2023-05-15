import cartopy
import numpy as np
import geospace as gs
from osgeo import gdal, ogr
from skimage.morphology import (remove_small_holes,
                                remove_small_objects)


def land_mask(out_file='/vsimem/land.tif', exclude_glacier=True, greenland=[126]):
    if gdal.Open(out_file) is not None:
        return out_file

    x_size, y_size, n_band, data_type = 3601, 1801, 1, gdal.GDT_Byte
    trans, srs = (-180.05, 0.1, 0.0, 90.05, 0.0, -0.1), 'EPSG:4326'

    # mask for lake
    f = '/vsimem/_lake.tif'
    shp = gs.shp_filter(cartopy.io.shapereader.natural_earth(name='lakes'),
                        filter_sql="scalerank = '0'")
    ds_shp = ogr.Open(shp)
    layer = ds_shp.GetLayer()
    gs.zeros_tif(f, x_size, y_size, n_band, data_type, trans, srs)
    ds = gdal.Open(f, gdal.GA_Update)
    gdal.RasterizeLayer(ds, [1], layer, burn_values=[1])
    ds = None
    lake = gdal.Open(f).ReadAsArray().astype(bool)

    # mask for Antarctica and Greenland
    f = '/vsimem/_glacier.tif'
    if exclude_glacier:
        sql = "','".join(str(i) for i in list(range(0, 8)) + greenland)
        shp = gs.shp_filter(cartopy.io.shapereader.natural_earth(name='land'),
                            filter_sql=f"FID IN ('{sql}')")
        ds_shp = ogr.Open(shp)
        layer = ds_shp.GetLayer()
        gs.zeros_tif(f, x_size, y_size, n_band, data_type, trans, srs)
        ds = gdal.Open(f, gdal.GA_Update)
        gdal.RasterizeLayer(ds, [1], layer, burn_values=[1],
                            options=['ALL_TOUCHED=TRUE'])
        ds = None
        glacier = gdal.Open(f).ReadAsArray().astype(bool)
    else:
        glacier = np.full(lake.shape, False)

    # mask for natureal earth land
    ds_shp = ogr.Open(cartopy.io.shapereader.natural_earth(name='land'))
    layer = ds_shp.GetLayer()
    gs.zeros_tif(out_file, x_size, y_size, n_band, data_type, trans, srs)
    ds = gdal.Open(out_file, gdal.GA_Update)
    gdal.RasterizeLayer(ds, [1], layer, burn_values=[1],
                        options=['ALL_TOUCHED=TRUE'])
    ds = None
    land = gdal.Open(out_file).ReadAsArray().astype(bool)

    # final mask
    ds = gdal.Open(out_file, gdal.GA_Update)
    valid = ~lake & ~glacier & land
    # valid = erosion(dilation(erosion(valid, disk(2)), disk(4)), disk(2))
    valid = remove_small_objects(valid, 30)
    valid = remove_small_holes(valid, 30)
    ds.WriteArray(valid)
    ds = None
    return out_file
