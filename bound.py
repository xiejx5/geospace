import numpy as np
from osgeo import gdal, ogr, osr
from .transform import geo2imagexy, ds_name


__all__ = ['bound_raster', 'bound_layers']


def proj_bound(ds, bound, bound_srs):
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

    # output SpatialReference
    if isinstance(bound_srs, osr.SpatialReference):
        outSpatialRef = bound_srs
    else:
        outSpatialRef = osr.SpatialReference()
        outSpatialRef.ImportFromProj4(bound_srs)

    # raster spatial reference
    ras_srs = osr.SpatialReference()
    ras_srs.ImportFromWkt(ds.GetProjection())

    # create the CoordinateTransformation
    trans = osr.CoordinateTransformation(outSpatialRef, ras_srs)
    trans_reverse = osr.CoordinateTransformation(ras_srs, outSpatialRef)

    # create a geometry from coordinates
    t = ds.GetGeoTransform()
    point = ogr.Geometry(ogr.wkbPoint)
    point_dx = ogr.Geometry(ogr.wkbPoint)
    point_dy = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(t[0] + ds.RasterXSize // 2 * t[1] + ds.RasterYSize // 2 *
                   t[2], t[3] + ds.RasterXSize // 2 * t[4] + ds.RasterYSize // 2 * t[5])
    point_dx.AddPoint(t[0] + (ds.RasterXSize // 2 + 1) * t[1] + ds.RasterYSize // 2 *
                      t[2], t[3] + (ds.RasterXSize // 2 + 1) * t[4] + ds.RasterYSize // 2 * t[5])
    point_dy.AddPoint(t[0] + ds.RasterXSize // 2 * t[1] + (ds.RasterYSize // 2 + 1) *
                      t[2], t[3] + ds.RasterXSize // 2 * t[4] + (ds.RasterYSize // 2 + 1) * t[5])
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
    bound = [win[:, 0:2].min(), win[:, 2:].min(),
             win[:, 0:2].max(), win[:, 2:].max()]
    return bound


def bound_raster(ds, bound, bound_srs="+proj=longlat +datum=WGS84 +ellps=WGS84"):
    ds, ras = ds_name(ds)
    t = ds.GetGeoTransform()
    x_min, y_min, x_max, y_max = proj_bound(ds, bound, bound_srs)
    ulX, ulY = geo2imagexy(ds, x_min, y_min)
    lrX, lrY = geo2imagexy(ds, x_max, y_max)
    clip_range = [min(ulX, lrX), min(ulY, lrY),
                  abs(ulX - lrX) + 1, abs(ulY - lrY) + 1]
    ul_lon = t[0] + t[1] * clip_range[0] + t[2] * clip_range[1]
    ul_lat = t[3] + t[4] * clip_range[0] + t[5] * clip_range[1]
    lr_lon = t[0] + t[1] * (clip_range[0] + clip_range[2]) + \
        t[2] * (clip_range[1] + clip_range[3])
    lr_lat = t[3] + t[4] * (clip_range[0] + clip_range[2]) + \
        t[5] * (clip_range[1] + clip_range[3])
    bound = [min(ul_lon, lr_lon), min(ul_lat, lr_lat),
             max(ul_lon, lr_lon), max(ul_lat, lr_lat)]

    # raster spatial reference as bound srs
    bound_srs = osr.SpatialReference()
    bound_srs.ImportFromWkt(ds.GetProjection())

    return bound, bound_srs.ExportToProj4()


def bound_layers(layers):
    # clip extent
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if not isinstance(layers, list):
        layers = [layers]
    if isinstance(layers[0], str):
        shp_in = [driver.Open(l) for l in layers]
        layers = [s.GetLayer() for s in shp_in]
    win = np.array([l.GetExtent() for l in layers])
    bound = [win[:, 0:2].min(), win[:, 2:].min(),
             win[:, 0:2].max(), win[:, 2:].max()]
    bound_srs = layers[0].GetSpatialRef()

    return bound, bound_srs.ExportToProj4()
