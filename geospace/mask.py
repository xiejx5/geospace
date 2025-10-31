import numpy as np
from osgeo import gdal, ogr
from geospace.utils import zeros_tif
from geospace.shape import shp_filter


def rounder(x, n=2):
    """Rounds a number to a specified number of significant figures.

    Args:
        x (float or int): The number to round.
        n (int, optional): The number of significant figures. Defaults to 2.

    Returns:
        float: The rounded number.
    """
    import math

    try:
        return round(x, -int(math.floor(math.log10(abs(math.modf(x)[0])))) + n - 1)
    except ValueError:
        return x


def shape_to_trans(y_size, x_size):
    """Calculates a GDAL geotransform from the shape of a global raster.

    Args:
        y_size (int): The number of rows (height).
        x_size (int): The number of columns (width).

    Returns:
        tuple: The GDAL geotransform.
    """
    # sizes=3600x1801 special case for ERA5
    res = rounder(360 / x_size)
    res = 360 / x_size if res < 0.01 else res
    trans = (
        -rounder(res * x_size / 2, 3),
        res,
        0.0,
        rounder(res * y_size / 2, 3),
        0.0,
        -res,
    )
    return trans


def land_mask(
    out_file='/vsimem/land.tif',
    sizes='3600x1800',
    greenland=[126],
    exclude_glacier=True,
):
    """Creates a land mask from Natural Earth data.

    This function generates a raster mask that separates land from water,
    with options to exclude glaciers and fill small holes.

    Args:
        out_file (str, optional): The path to the output mask file.
                                  Defaults to '/vsimem/land.tif'.
        sizes (str, optional): The dimensions of the output mask (e.g., '3600x1800').
                               Defaults to '3600x1800'.
        greenland (list, optional): A list of FIDs for Greenland polygons to exclude.
                                    Defaults to [126].
        exclude_glacier (bool, optional): If True, exclude glaciers. Defaults to True.

    Returns:
        str: The path to the output mask file.
    """
    import re
    import math
    from cartopy.io.shapereader import natural_earth
    from skimage.morphology import remove_small_holes, remove_small_objects

    if gdal.Open(out_file) is not None:
        return out_file

    x_size, y_size = (int(size) for size in re.findall(r'\d+', sizes))
    trans = shape_to_trans(y_size, x_size)
    n_band, data_type, srs, res = 1, gdal.GDT_Byte, 'EPSG:4326', trans[1]

    # mask for lake
    f = '/vsimem/_lake.tif'
    ds_shp = shp_filter(natural_earth(name='lakes'), "scalerank = '0'")
    layer = ds_shp.GetLayer()
    zeros_tif(f, x_size, y_size, n_band, data_type, trans, srs)
    ds = gdal.Open(f, gdal.GA_Update)
    gdal.RasterizeLayer(ds, [1], layer, burn_values=[1])
    ds = None
    lake = gdal.Open(f).ReadAsArray().astype(bool)

    # mask for Antarctica and Greenland
    f = '/vsimem/_glacier.tif'
    if exclude_glacier:
        sql = "','".join(str(i) for i in list(range(0, 8)) + greenland)
        ds_shp = shp_filter(natural_earth(name='land'), f"FID IN ('{sql}')")
        layer = ds_shp.GetLayer()
        zeros_tif(f, x_size, y_size, n_band, data_type, trans, srs)
        ds = gdal.Open(f, gdal.GA_Update)
        gdal.RasterizeLayer(
            ds, [1], layer, burn_values=[1], options=['ALL_TOUCHED=TRUE']
        )
        ds = None
        glacier = gdal.Open(f).ReadAsArray().astype(bool)
    else:
        glacier = np.full(lake.shape, False)

    # mask for natureal earth land
    ds_shp = ogr.Open(natural_earth(name='land'))
    layer = ds_shp.GetLayer()
    zeros_tif(out_file, x_size, y_size, n_band, data_type, trans, srs)
    ds = gdal.Open(out_file, gdal.GA_Update)
    gdal.RasterizeLayer(ds, [1], layer, burn_values=[1], options=['ALL_TOUCHED=TRUE'])
    ds = None
    land = gdal.Open(out_file).ReadAsArray().astype(bool)

    # final mask
    ds = gdal.Open(out_file, gdal.GA_Update)
    valid = ~lake & ~glacier & land
    # valid = erosion(dilation(erosion(valid, disk(2)), disk(4)), disk(2))
    # fill or remove 30 grids for 0.1
    thres = math.ceil(0.3 / (res**2))

    # morphology
    valid = remove_small_objects(valid, thres)
    valid = remove_small_holes(valid, thres)

    # save to tif
    ds.WriteArray(valid)
    ds = None
    return out_file
