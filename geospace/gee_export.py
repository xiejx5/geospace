import os
import ee
import zipfile
import requests
import pandas as pd


# gee initialization
def gee_initial():
    try:
        ee.Initialize()
    except Exception:
        ee.Authenticate()
        ee.Initialize()


def gee_export_tif(image, filename, crs=None, crs_transform=None, scale=None, region=None, file_per_band=False):
    """Export as tif, must be the original image, instead of the reprojected one

    Args:
        fc (ee.FeatureCollection): the spatial scope of the exported tif
        image (ee.Image): the image that would be exported
        filename (string): exported path
        crs (str, optional): A default CRS string to use for any bands that do not explicitly specify one
        crs_transform ([type]): control the spatial resolution and alignment
        region (object, optional): A polygon specifying a region to download
        file_per_band (bool, optional): Whether to produce a different GeoTIFF per band
    """
    if not isinstance(image, ee.Image):
        print("The image must be an ee.Image.")
        return

    filename = os.path.abspath(filename)
    basename = os.path.basename(filename)
    name = os.path.splitext(basename)[0]
    filetype = os.path.splitext(basename)[1][1:].lower()
    filename_zip = filename.replace(".tif", ".zip")

    if filetype != "tif":
        print("The filename must end with .tif")
        return

    try:
        print("Generating URL ...")
        params = {"name": name, "filePerBand": file_per_band}
        if region is None:
            region = image.geometry().getInfo()
        params["region"] = region
        if crs is not None:
            params["crs"] = crs
        if crs_transform is not None:
            params["crs_transform"] = crs_transform
        elif scale is not None:
            params["scale"] = scale
        else:
            params["scale"] = image.projection().nominalScale()

        url = image.getDownloadURL(params)
        print(f"Downloading data from {url}\nPlease wait ...")
        r = requests.get(url, stream=True)

        if r.status_code != 200:
            print("An error occurred while downloading.")
            return

        with open(filename_zip, "wb") as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)

    except Exception as e:
        print("An error occurred while downloading.")
        print(e)
        return

    try:
        with zipfile.ZipFile(filename_zip) as z:
            z.extractall(os.path.dirname(filename))
        os.remove(filename_zip)

        if file_per_band:
            print(f"Data downloaded to {os.path.dirname(filename)}")
        else:
            print(f"Data downloaded to {filename}")
    except Exception as e:
        print(e)

    # ee.batch.Export.image.toDrive(
    #     **{
    #         'image': image,
    #         'description': asset,
    #         'fileNamePrefix': asset,  # this is the name actually
    #         'folder': 'SoilGrids',
    #         'region': fc.geometry(),
    #         'crs': 'EPSG:4326',
    #         'crsTransform': crs_transform,
    #         'maxPixels': 1e13
    #     }).start()

    # ee.batch.Export.image.toDrive(
    #     **{
    #         'image': image,
    #         'description': asset,
    #         'fileNamePrefix': asset,  # this is the name actually
    #         'folder': 'SoilGrids',
    #         'region': fc.geometry(),
    #         'crs': 'EPSG:4326',
    #         'crsTransform': crs_transform,
    #         'maxPixels': 1e13
    #     }).start()


def gee_export_csv(fc, image, scale_enlarge=1, fields=['ORDER', '.*mean'], return_url=False):
    """export a csv containing the basin average value

    Args:
        fc (ee.FeatureCollection): e.g. basins
        image (ee.Image): e.g. DEM

    Returns:
        DataFrame: it has fields of 'ORDER' and 'mean'
    """
    # export as csv
    scale = image.projection().nominalScale().multiply(scale_enlarge)

    if image.projection().crs().getInfo() != 'EPSG:4326':
        image = image.reproject('EPSG:4326', None, scale)

    means = image.reduceRegions(**{
        'collection': fc,
        'reducer': ee.Reducer.mean(),
        'scale': scale,
    })
    url = means.select(fields, retainGeometry=False).getDownloadURL(filetype='csv')
    if return_url:
        return url
    else:
        return pd.read_csv(url)
