import numpy as np
from osgeo import ogr
from geospace._const import WGS84
from geospace.utils import geo2imagexy, ds_name
from geospace.projection import read_srs, coord_trans


def _prj_bound(ds, bound, bound_srs):
    x_min, y_min, x_max, y_max = bound
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(min(x_min, x_max), min(y_min, y_max))
    ring.AddPoint(max(x_min, x_max), min(y_min, y_max))
    ring.AddPoint(max(x_min, x_max), max(y_min, y_max))
    ring.AddPoint(min(x_min, x_max), max(y_min, y_max))
    ring.AddPoint(min(x_min, x_max), min(y_min, y_max))

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    trans = coord_trans(bound_srs, ds)
    trans_reverse = coord_trans(ds, bound_srs)

    # create a geometry from coordinates
    t = ds.GetGeoTransform()
    point = ogr.Geometry(ogr.wkbPoint)
    point_dx = ogr.Geometry(ogr.wkbPoint)
    point_dy = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(
        t[0] + ds.RasterXSize // 2 * t[1] + ds.RasterYSize // 2 * t[2],
        t[3] + ds.RasterXSize // 2 * t[4] + ds.RasterYSize // 2 * t[5],
    )
    point_dx.AddPoint(
        t[0] + (ds.RasterXSize // 2 + 1) * t[1] + ds.RasterYSize // 2 * t[2],
        t[3] + (ds.RasterXSize // 2 + 1) * t[4] + ds.RasterYSize // 2 * t[5],
    )
    point_dy.AddPoint(
        t[0] + ds.RasterXSize // 2 * t[1] + (ds.RasterYSize // 2 + 1) * t[2],
        t[3] + ds.RasterXSize // 2 * t[4] + (ds.RasterYSize // 2 + 1) * t[5],
    )
    point.Transform(trans_reverse)
    point_dx.Transform(trans_reverse)
    point_dy.Transform(trans_reverse)
    dx = abs(point.GetPoint()[0] - point_dx.GetPoint()[0])
    dy = abs(point.GetPoint()[1] - point_dy.GetPoint()[1])

    # density geom
    geom = ogr.CreateGeometryFromWkb(poly.ExportToWkb())
    geom.Segmentize(min(dx, dy) / 2)
    geom.Transform(trans)

    # get boundary
    win = np.array([geom.GetEnvelope()])
    bound = [win[:, 0:2].min(), win[:, 2:].min(), win[:, 0:2].max(), win[:, 2:].max()]
    return bound


def _enlarge_bound(ds, x_min, y_min, x_max, y_max):
    t = ds.GetGeoTransform()
    ulX, ulY = geo2imagexy(ds, x_min, y_min)
    lrX, lrY = geo2imagexy(ds, x_max, y_max)
    clip_range = [
        np.min([ulX, lrX], axis=0),
        np.min([ulY, lrY], axis=0),
        (np.abs(ulX - lrX) + 1),
        (np.abs(ulY - lrY) + 1),
    ]
    ul_lon = t[0] + t[1] * clip_range[0] + t[2] * clip_range[1]
    ul_lat = t[3] + t[4] * clip_range[0] + t[5] * clip_range[1]
    lr_lon = (
        t[0]
        + t[1] * (clip_range[0] + clip_range[2])
        + t[2] * (clip_range[1] + clip_range[3])
    )
    lr_lat = (
        t[3]
        + t[4] * (clip_range[0] + clip_range[2])
        + t[5] * (clip_range[1] + clip_range[3])
    )
    bound = [
        np.min([ul_lon, lr_lon], axis=0),
        np.min([ul_lat, lr_lat], axis=0),
        np.max([ul_lon, lr_lon], axis=0),
        np.max([ul_lat, lr_lat], axis=0),
    ]

    return bound, clip_range


def bound_raster(ds, bound, bound_srs=WGS84):
    ds, _ = ds_name(ds)
    x_min, y_min, x_max, y_max = _prj_bound(ds, bound, bound_srs)
    bound = _enlarge_bound(ds, x_min, y_min, x_max, y_max)[0]

    return bound, read_srs(ds).ExportToProj4()


def bound_layers(layers):
    # clip extent
    if not isinstance(layers, list):
        layers = [layers]
    if isinstance(layers[0], str):
        shp_in = [ogr.Open(lyr) for lyr in layers]
        layers = [s.GetLayer() for s in shp_in]
    win = np.array([lyr.GetExtent() for lyr in layers])
    bound = [win[:, 0:2].min(), win[:, 2:].min(), win[:, 0:2].max(), win[:, 2:].max()]

    return bound, read_srs(layers).ExportToProj4()


def grid_bound(ds, regions):
    """return gridded bound whose resolution exactly matches that of ds

    Args:
        ds (gdal.dataset): dataset or path
        regions (list or str): regions of interest

    Returns:
        bound: spatial range of the regions in the ds
        crs_transform: like bound but used in GEE
    """
    ds = ds_name(ds)[0]
    t = ds.GetGeoTransform()
    bound, bound_srs = bound_layers(regions)
    bound = bound_raster(ds, bound, bound_srs)[0]
    crs_transform = [t[1], 0, bound[0], 0, t[5], bound[-1]]
    return bound, crs_transform
