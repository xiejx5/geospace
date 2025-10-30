import os
import numpy as np
from pathlib import Path
from osgeo import gdal, ogr
from functools import partial
from geospace.projection import read_srs
from geospace.boundary import _enlarge_bound
from geospace.spatial_calc import area_per_row
from geospace._const import WGS84, CREATION, N_CPU
from geospace.shape import shp_projection, shp_filter
from geospace.utils import zeros_tif, block_write, context_file, ds_name, rep_name


# Now Read the large raster block by block, expecting to see no increase of
# memory as we loop through the blocks.
def _map_burn(arrs, burn_data):
    mask = np.broadcast_to(np.logical_not(burn_data.astype(bool)), arrs[0].shape)
    rect_data = np.ma.masked_where(mask, arrs[0])
    return rect_data.filled()


def _in_shape(rect_trans, enlarge, n_x, n_y, outLayer, rasterize_option):
    # set geotransform for enlarged tif
    poly_trans = list(rect_trans)
    poly_trans[1] = poly_trans[1] / enlarge
    poly_trans[5] = poly_trans[5] / enlarge

    # Rasterize
    poly_file = '/vsimem/_poly.tif'
    zeros_tif(
        poly_file,
        n_x * enlarge,
        n_y * enlarge,
        1,
        gdal.GDT_Byte,
        poly_trans,
        outLayer.GetSpatialRef(),
        nodata=2,
    )
    poly_ds = gdal.Open(poly_file, gdal.GA_Update)
    gdal.RasterizeLayer(
        poly_ds, [1], outLayer, burn_values=[1], options=rasterize_option
    )
    poly_data = poly_ds.ReadAsArray()

    # https://towardsdatascience.com/efficiently-splitting-an-image-into-tiles-in-python-using-numpy-d1bf0dd7b6f7
    if enlarge == 1:
        count = poly_data
    else:
        count = poly_data.reshape(n_y, enlarge, n_x, enlarge).swapaxes(1, 2)
        count = count.sum(axis=(-2, -1), dtype=np.uint16)
    return count


def _clip(
    ds,
    outLayer,
    out_file,
    ext='',
    enlarge=None,
    save_cache=False,
    reuse_cache=False,
    rasterize_option=['ALL_TOUCHED=TRUE'],
):
    # auto downscale the raster if resolution > 50 meter
    t = ds.GetGeoTransform()
    if enlarge is None:
        if 'PROJCS' in ds.GetProjection():
            enlarge = 1 if abs(t[1]) < 50 else 10
        else:
            enlarge = 1 if abs(t[1]) < 0.0005 else 10

    # get no data and Spatial Reference
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    srs = outLayer.GetSpatialRef()

    # get shp extent in form of raster grids
    x_min, x_max, y_min, y_max = outLayer.GetExtent()
    bound, clip_range = _enlarge_bound(ds, x_min, y_min, x_max, y_max)
    n_x, n_y = int(clip_range[2]), int(clip_range[3])

    if out_file is not None:
        # clip with rectangle out_file, use Warp instead of Translate
        # to project longitude from (0, 360) to (-180, 180)
        option = gdal.WarpOptions(
            multithread=True,
            outputBounds=bound,
            srcSRS=srs,
            dstSRS=srs,
            creationOptions=CREATION,
            dstNodata=nodata,
            xRes=t[1],
            yRes=t[5],
            srcNodata=nodata,
            resampleAlg=gdal.GRA_NearestNeighbour,
        )
        ds_rect = gdal.Warp(out_file, ds, options=option)
        # mask value outside the shapefile
        is_in = _in_shape(
            ds_rect.GetGeoTransform(), 1, n_x, n_y, outLayer, rasterize_option
        )
        block_write(ds_rect, _map_burn, is_in)
        return ds_rect, is_in

    # get a rectangle containing the shapefile
    rect_trans = list(t)
    rect_trans[0], rect_trans[3] = bound[0], bound[3]
    # longitude from (0, 360) to (-180, 180)
    xoff = clip_range[0] if clip_range[0] > 0 else clip_range[0] + ds.RasterXSize
    xoff, yoff = int(xoff), int(clip_range[1])
    # shape cross the prime meridian
    if xoff + n_x > ds.RasterXSize:
        rect = np.concatenate(
            [
                ds.ReadAsArray(xoff, yoff, ds.RasterXSize - xoff, n_y),
                ds.ReadAsArray(0, yoff, xoff + n_x - ds.RasterXSize, n_y),
            ],
            axis=-1,
        )
    else:
        rect = ds.ReadAsArray(xoff, yoff, n_x, n_y)

    # paths of cached _burn.tif
    if save_cache:
        reuse_cache = True
        cache_dir = 'cache'
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
        burn_file = os.path.join(cache_dir, str(ext) + '_burn.tif')
    else:
        burn_file = os.path.join('/vsimem/', str(ext) + '_burn.tif')

    # calculate shapefile intersection area in each grid
    burn_ds = gdal.Open(burn_file)
    if (
        reuse_cache
        and (burn_ds is not None)
        and np.array_equal(rect_trans, burn_ds.GetGeoTransform())
        and srs.ExportToProj4() == burn_ds.GetSpatialRef().ExportToProj4()
        and n_x == burn_ds.RasterXSize
        and n_y == burn_ds.RasterYSize
    ):
        burn_data = burn_ds.ReadAsArray()
    else:
        # get counts of shape intersection for each grid
        count = _in_shape(rect_trans, enlarge, n_x, n_y, outLayer, rasterize_option)

        # calculate area for each grid
        row_area = area_per_row(rect_trans, srs.ExportToWkt(), n_y)
        burn_data = count * row_area.reshape([n_y, -1]) / enlarge / enlarge

        # save the burn_data into tif
        if reuse_cache:
            burn_ds = gdal.GetDriverByName('GTiff').Create(
                burn_file, n_x, n_y, 1, gdal.GDT_Float64, CREATION
            )
            burn_ds.SetGeoTransform(tuple(rect_trans))
            burn_ds.SetSpatialRef(srs)
            burn_ds.GetRasterBand(1).SetNoDataValue(0)
            burn_ds.WriteArray(burn_data)
            burn_ds = None

    return rect, burn_data


def extract(ras, shp, out_path=None, ras_srs=WGS84, nodata=None, **kwargs):
    ds, ras = ds_name(ras)
    if out_path is None:
        out_file = None
    else:
        out_file = context_file(ds.GetDescription(), out_path)
        if os.path.exists(out_file):
            return out_file

    # set projection
    inDataset = shp if isinstance(shp, gdal.Dataset) else ogr.Open(shp)
    inLayer = inDataset.GetLayer()
    if (
        read_srs([ds, ras_srs]).ExportToProj4()
        == inLayer.GetSpatialRef().ExportToProj4()
    ):
        outLayer = inLayer
    else:
        outDataSet = shp_projection(inLayer, out_srs=read_srs([ds, ras_srs]))
        outLayer = outDataSet.GetLayer()

    # set no data
    if ds.GetRasterBand(1).GetNoDataValue() is not None:
        nodata = ds.GetRasterBand(1).GetNoDataValue()
    elif nodata is not None:
        ds.GetRasterBand(1).SetNoDataValue(nodata)
    else:
        raise (ValueError('nodata must be initialized'))

    # clip with the whole shapefile
    rect, burn_data = _clip(ds, outLayer, out_file, **kwargs)
    if out_file is not None:
        return rect.GetDescription()

    # area weighted statistics
    not_in = np.broadcast_to(np.logical_not(burn_data.astype(bool)), rect.shape)
    mask = (rect == nodata) | not_in | (~np.isfinite(rect))
    arr = np.ma.masked_array(rect, mask)
    arr = arr[np.newaxis, :, :] if arr.ndim != 3 else arr

    return np.ma.average(
        arr.reshape(arr.shape[0], -1), weights=burn_data.ravel(), axis=1
    ).filled(np.nan)


def basin_average_worker(rasters, shp, is_unique, s, t, field, sel, **kwargs):
    where = f'{field}={f"'{sel.replace("'", "''")}'" if isinstance(sel, str) else sel}'
    ds_shp = shp_filter(shp, where)
    one_out = np.full(t[-1], np.nan)
    for i, ras in enumerate(rasters):
        kwargs['reuse_cache'] = False if is_unique[i] else True
        one_out[s[i] : t[i]] = extract(ras, ds_shp, ext=sel, **kwargs)
    return one_out


def basin_average(rasters, shp, field='STAID', sel=None, **kwargs):
    import tqdm
    import pandas as pd
    from multiprocessing import get_context

    if isinstance(rasters, (str, Path)):
        rasters = [str(rasters)]
    else:
        rasters = [str(ras) for ras in rasters]

    if sel is None:
        ds = ogr.Open(shp)
        layer = ds.GetLayer()
        layer.ResetReading()
        sel = [feature.GetFID() for feature in layer]
        field = 'FID'
    if isinstance(sel, (str, int)):
        sel = [sel]

    arr = np.array([gdal.Open(ras).GetGeoTransform() for ras in rasters])
    sort_idxs = np.lexsort((arr[:, 0], arr[:, 3], arr[:, 1], arr[:, 5]))
    sort_rasters = np.array(rasters)[sort_idxs]
    _, idxs, counts = np.unique(
        arr[:, [1, 5]][sort_idxs], axis=0, return_index=True, return_counts=True
    )
    is_unique = np.full(len(rasters), False)
    is_unique[idxs[counts == 1]] = True

    sort_names, s, t, inverse = rep_name(sort_rasters, sort_idxs=sort_idxs)

    cpu_used = max(min(N_CPU, len(sel)), 1)
    if cpu_used == 1 or kwargs.pop('parallel', True) == False:
        output = list(
            tqdm.tqdm(
                (
                    partial(
                        basin_average_worker,
                        sort_rasters,
                        shp,
                        is_unique,
                        s,
                        t,
                        field,
                        **kwargs,
                    )(f)
                    for f in sel
                ),
                total=len(sel),
            )
        )
    else:
        with get_context('spawn').Pool(cpu_used) as p:
            output = list(
                tqdm.tqdm(
                    p.imap(
                        partial(
                            basin_average_worker,
                            sort_rasters,
                            shp,
                            is_unique,
                            s,
                            t,
                            field,
                            **kwargs,
                        ),
                        sel,
                    ),
                    total=len(sel),
                )
            )

    # use column to store array for performance
    # df = pd.DataFrame(output, columns=sort_names, index=sel)
    df = pd.DataFrame(np.vstack(output).T, columns=sel, index=sort_names)

    return df.iloc[inverse, :]
