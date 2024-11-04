import re
import math
import cartopy
import numpy as np
import geospace as gs
from osgeo import gdal, ogr
from skimage.morphology import (remove_small_holes,
                                remove_small_objects)


def rounder(x, n=2):
    """Round a value to n significant figures.

    If x is an integer or have only zeros after the decimal point, no need to round.
    If x is a float, keep only the first n digits other than zero after the decimal point
    """
    try:
        return round(x, -int(math.floor(math.log10(abs(math.modf(x)[0])))) + n - 1)
    except ValueError:
        return x


def shape_to_trans(y_size, x_size):
    # sizes=3600x1801 special case for ERA5
    res = rounder(360 / x_size)
    res = 360 / x_size if res < 0.01 else res
    trans = (-rounder(res * x_size / 2, 3), res, 0.0,
             rounder(res * y_size / 2, 3), 0.0, -res)
    return trans


def land_mask(out_file='/vsimem/land.tif', sizes='3600x1800', greenland=[126], exclude_glacier=True):
    if gdal.Open(out_file) is not None:
        return out_file

    x_size, y_size = (int(size) for size in re.findall(r'\d+', sizes))
    trans = shape_to_trans(y_size, x_size)
    n_band, data_type, srs, res = 1, gdal.GDT_Byte, 'EPSG:4326', trans[1]

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
    # fill or remove 30 grids for 0.1
    thres = math.ceil(0.3 / (res**2))
    valid = remove_small_objects(valid, thres)
    valid = remove_small_holes(valid, thres)
    ds.WriteArray(valid)
    ds = None
    return out_file
