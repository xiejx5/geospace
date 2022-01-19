from pathlib import Path
from typing import List
from osgeo import gdal, ogr, osr


def read_srs(srs):
    if isinstance(srs, osr.SpatialReference):
        return srs
    elif isinstance(srs, str):
        if Path(srs).is_file() or 'vsimem' in str(Path(srs).parent):
            ds = gdal.Open(srs)
            if ds is not None:
                return osr.SpatialReference(wkt=ds.GetProjection())
            ds = ogr.Open(srs)
            if ds is not None:
                layer = ds.GetLayer()
                return layer.GetSpatialRef()
        else:
            SpatialRef = osr.SpatialReference(wkt=None)
            if SpatialRef.ImportFromProj4(srs) == 0:
                return SpatialRef
            if SpatialRef.ImportFromWkt(srs) == 0:
                return SpatialRef
    elif isinstance(srs, gdal.Dataset) or isinstance(srs, ogr.DataSource):
        return read_srs(srs.GetDescription())
    elif isinstance(srs, ogr.Layer):
        return srs.GetSpatialRef()
    elif isinstance(srs, List):
        for s in srs:
            SpatialRef = read_srs(s)
            if SpatialRef is not None and SpatialRef.ExportToWkt() != '':
                return SpatialRef
        raise(ValueError("srs must be set"))
    return None


def coord_trans(in_srs, out_srs):
    inSpatialRef, outSpatialRef = read_srs(in_srs), read_srs(out_srs)
    # create the CoordinateTransformation
    if int(gdal.__version__[0]) >= 3:
        outSpatialRef.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    return coordTrans
