import numpy as np
from collections.abc import Iterable


__all__ = ['block_write']


def block_write(in_ds, in_bands, out_band, map_fun,
                map_args=None, map_kwargs=None):
    if isinstance(in_bands, Iterable):
        block_write_multi(in_ds, in_bands, out_band, map_fun,
                          map_args=map_args, map_kwargs=map_kwargs)
    else:
        block_write_one(in_ds, in_bands, out_band, map_fun,
                        map_args=map_args, map_kwargs=map_kwargs)


def block_write_one(in_ds, in_band, out_band, map_fun,
                    map_args=None, map_kwargs=None):
    block_xsize, block_ysize = in_band.GetBlockSize()
    for b_y, yoff in enumerate(range(0, in_ds.RasterYSize, block_ysize)):
        for b_x, xoff in enumerate(range(0, in_ds.RasterXSize, block_xsize)):
            win_xsize, win_ysize = in_band.GetActualBlockSize(b_x, b_y)

            in_data = in_band.ReadAsArray(
                xoff=xoff, yoff=yoff, win_xsize=win_xsize, win_ysize=win_ysize)
            in_data = np.ma.masked_array(in_data,
                                         mask=in_data == in_band.GetNoDataValue(),
                                         fill_value=out_band.GetNoDataValue())

            if map_args and map_kwargs:
                out_data = map_fun(in_data, *map_args, **map_kwargs)
            elif map_args:
                out_data = map_fun(in_data, *map_args)
            elif map_kwargs:
                out_data = map_fun(in_data, **map_kwargs)
            else:
                out_data = map_fun(in_data)

            out_band.WriteArray(out_data, xoff=xoff, yoff=yoff)


def block_write_multi(in_ds, in_bands, out_band, map_fun,
                      map_args=None, map_kwargs=None):
    first_band = in_bands[0]
    block_xsize, block_ysize = first_band.GetBlockSize()
    for b_y, yoff in enumerate(range(0, in_ds.RasterYSize, block_ysize)):
        for b_x, xoff in enumerate(range(0, in_ds.RasterXSize, block_xsize)):
            win_xsize, win_ysize = first_band.GetActualBlockSize(b_x, b_y)
            no_data = out_band.GetNoDataValue()

            in_datas = []
            for in_band in in_bands:
                in_data = in_band.ReadAsArray(
                    xoff=xoff, yoff=yoff, win_xsize=win_xsize, win_ysize=win_ysize)
                in_data = np.ma.masked_array(in_data,
                                             mask=in_data == in_band.GetNoDataValue(),
                                             fill_value=no_data)
                in_datas.append(in_data)

            if map_args and map_kwargs:
                out_data = map_fun(in_datas, *map_args, **map_kwargs)
            elif map_args:
                out_data = map_fun(in_datas, *map_args)
            elif map_kwargs:
                out_data = map_fun(in_datas, **map_kwargs)
            else:
                out_data = map_fun(in_datas)

            out_band.WriteArray(out_data, xoff=xoff, yoff=yoff)
