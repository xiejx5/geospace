import os
import numpy as np
from osgeo import gdal
from pathlib import Path
from geospace.raster import mosaic
from collections.abc import Iterable
from geospace._const import CREATION, N_CPU
from geospace.utils import ds_name, context_file


def band_map(i, ras_multi, band_idx_multi, calc_arg, out_file):
    from osgeo_utils.gdal_calc import Calc

    ds_ras = gdal.Open(ras_multi[0])
    band_name = ds_ras.GetRasterBand(int(band_idx_multi[0])).GetDescription()
    tem_name = band_name if band_name else f'_temp_{Path(out_file).stem}_{str(i)}'
    tem_file = f'{Path(out_file).parent}/{tem_name}.tif'

    if os.path.exists(tem_file):
        return tem_file

    ras_args = {chr(i + 65): ras for i, ras in enumerate(ras_multi)}
    band_args = {
        chr(i + 65) + '_band': int(band_idx)
        for i, band_idx in enumerate(band_idx_multi)
    }
    input_args = {**ras_args, **band_args}

    Calc(calc_arg, tem_file, creation_options=CREATION, quiet=True, **input_args)

    return tem_file


def check_iter(ds_multi, calc_args, band_idxs):
    iter_ds_multi = isinstance(ds_multi, Iterable) and not isinstance(ds_multi, str)
    iter_calc_args = isinstance(calc_args, Iterable) and not isinstance(calc_args, str)
    if band_idxs is not None:
        if iter_ds_multi:
            iter_band_idxs = isinstance(band_idxs[0], Iterable)
        else:
            iter_band_idxs = isinstance(band_idxs, Iterable)
    else:
        iter_band_idxs = False

    return iter_ds_multi, iter_calc_args, iter_band_idxs


def broadcast_args(ds_multi, calc_args, band_idxs):
    iter_ds_multi, iter_calc_args, iter_band_idxs = check_iter(
        ds_multi, calc_args, band_idxs
    )

    if iter_ds_multi:
        ds = ds_multi[0]
    else:
        ds = ds_multi

    ds, ras = ds_name(ds)
    if band_idxs is not None:
        if iter_band_idxs and iter_calc_args:
            if len(band_idxs) != len(calc_args):
                raise Exception('length of band list not equal to that of calc args')
        elif iter_band_idxs:
            calc_args = [calc_args] * len(band_idxs)
        elif iter_calc_args:
            band_idxs = [band_idxs] * len(calc_args)
        else:
            calc_args = [calc_args]
            band_idxs = [band_idxs]
    else:
        n_band = ds.RasterCount
        if iter_calc_args:
            if len(calc_args) != n_band:
                raise Exception('calc args length not equal to band counts')
        else:
            calc_args = [calc_args] * n_band

        if iter_ds_multi:
            band_idxs = np.repeat(
                np.arange(1, n_band + 1, dtype=int), len(ds_multi)
            ).reshape(-1, len(ds_multi))
        else:
            band_idxs = np.arange(1, n_band + 1, dtype=int)

    if not iter_ds_multi:
        ds_multi = [ras]
        band_idxs = np.array(band_idxs).reshape(len(band_idxs), 1)

    return ds_multi, calc_args, band_idxs


def map_calc(ds_multi, calc_args, out_path, band_idxs=None, parallel=True):
    import tqdm
    from multiprocessing import get_context

    if isinstance(ds_multi, Iterable) and not isinstance(ds_multi, str):
        ds_multi = [str(ds) for ds in ds_multi]
        ds = ds_multi[0]
    else:
        ds_multi = str(ds_multi)
        ds = ds_multi

    ds, ras = ds_name(ds)
    out_file = context_file(ras, out_path)

    if os.path.exists(out_file):
        return out_file

    ds_multi, calc_args, band_idxs = broadcast_args(ds_multi, calc_args, band_idxs)

    n = len(calc_args)
    args = zip(
        np.arange(1, n + 1, dtype=int),
        [ds_multi] * n,
        band_idxs,
        calc_args,
        [out_file] * n,
    )

    if parallel and (n > 1):
        with get_context('spawn').Pool(min(N_CPU, n)) as p:
            tem_files = p.starmap(band_map, tqdm.tqdm(args, total=n))
    else:
        tem_files = []
        for arg in args:
            tem_files.append(band_map(*arg))

    if len(tem_files) == 1:
        os.rename(tem_files[0], out_file)
    else:
        mosaic(tem_files, out_file, separate=True)
        [os.remove(f) for f in tem_files]

    return out_file
