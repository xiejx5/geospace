import os


def connected_to_internet(url='http://www.google.com/', timeout=1):
    import requests

    try:
        _ = requests.head(url, timeout=timeout)
        return True
    except requests.ConnectionError:
        print('Failed to access GEE')
    return False


# gee initialization
def gee_init(proxy='http://127.0.0.1:7890'):
    import ee

    if not connected_to_internet():
        print('Set proxy to 127.0.0.1:7890')
        os.environ['HTTP_PROXY'] = proxy
        os.environ['HTTPS_PROXY'] = proxy
    try:
        ee.Initialize()
    except Exception:
        ee.Authenticate()
        ee.Initialize()


def gee_export_tif(image, filename, timeout=300, **params):
    """Export as tif, must be the original image, instead of the reprojected one

    Args:
        fc (ee.FeatureCollection): the spatial scope of the exported tif
        image (ee.Image): the image that would be exported
        filename (string): exported path
        crs (str, optional): A default CRS string to use for any bands that do not explicitly specify one
        crs_transform ([type]): control the spatial resolution and alignment, e.g., [1, 0, -180, 0, -1, 90]
        region (object, optional): A polygon specifying a region to download
        filePerBand (bool, optional): Whether to produce a different GeoTIFF per band

        details are in ee.Image.getDownloadURL()
    """
    import ee
    import requests
    import zipfile

    if not isinstance(image, ee.Image):
        print('The image must be an ee.Image.')
        return

    filename = os.path.abspath(filename)
    basename = os.path.basename(filename)
    name = os.path.splitext(basename)[0]
    filetype = os.path.splitext(basename)[1][1:].lower()
    filename_zip = filename.replace('.tif', '.zip')

    if filetype != 'tif':
        print('The filename must end with .tif')
        return

    try:
        print('Generating URL ...')
        params['name'] = name
        params['filePerBand'] = params.pop('filePerBand', False)
        if params.get('region') == 'global':
            params['region'] = image.geometry()

        try:
            url = image.getDownloadURL(params)
        except Exception as e:
            print('An error occurred while downloading.')
            print(e)
            return
        print(f'Downloading data from {url}\nPlease wait ...')
        r = requests.get(url, stream=True, timeout=timeout)

        if r.status_code != 200:
            print('An error occurred while downloading.')
            return

        with open(filename_zip, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)

    except Exception as e:
        print('An error occurred while downloading.')
        print(r.json()['error']['message'])
        return

    try:
        with zipfile.ZipFile(filename_zip) as z:
            z.extractall(os.path.dirname(filename))
        os.remove(filename_zip)

        if params['filePerBand']:
            print(f'Data downloaded to {os.path.dirname(filename)}')
        else:
            print(f'Data downloaded to {filename}')
    except Exception as e:
        print(e)


def gee_to_drive(image, **params):
    """wrapper for GEE Export.image.toDrive

    Args:
        image: The image to be exported.
        description: Human-readable name of the task.
        folder: The name of a unique folder in your Drive account to
            export into. Defaults to the root of the drive.
        dimensions: The dimensions of the exported image. Takes either a
            single positive integer as the maximum dimension or "WIDTHxHEIGHT"
            where WIDTH and HEIGHT are each positive integers.
        region: The lon,lat coordinates for a LinearRing or Polygon
            specifying the region to export. Can be specified as a nested
            lists of numbers or a serialized string. Defaults to the image's
            region.
        crs: The coordinate reference system of the exported image's
            projection. Defaults to the image's default projection.
        crsTransform: A comma-separated string of 6 numbers describing
            the affine transform of the coordinate reference system of the
            exported image's projection, in the order: xScale, xShearing,
            xTranslation, yShearing, yScale and yTranslation. Defaults to
            the image's native CRS transform.
        maxPixels: The maximum allowed number of pixels in the exported
            image. The task will fail if the exported region covers more
            pixels in the specified projection. Defaults to 100,000,000.

        details are in Export.image.toDrive()
    """
    import ee

    params['image'] = image
    params['description'] = params.pop(
        'description', image.bandNames().get(0).getInfo()
    )
    params['crs'] = params.pop('crs', 'EPSG:4326')
    params['maxPixels'] = params.pop('maxPixels', 1e13)

    if params.get('crsTransform') == 'ERA5_LAND':
        params['crsTransform'] = [0.1, 0, -180.05, 0, -0.1, 90.05]

    if params.get('region') == 'global':
        params['region'] = ee.Geometry.Rectangle(
            [-180, -90, 180, 90], geodesic=False, proj='EPSG:4326'
        )
        if params.get('dimensions') == 'default':
            params['dimensions'] = '360x180'
    elif params.get('region') == 'nonATA':
        params['region'] = ee.Geometry.Rectangle(
            [-180, -60, 180, 89], geodesic=False, proj='EPSG:4326'
        )
        if params.get('dimensions') == 'default':
            params['dimensions'] = '360x149'

    ee.batch.Export.image.toDrive(**params).start()


def gee_export_csv(
    fc, image, fields=['STAID', '.*mean'], return_url=False, to_drive=None, **kwargs
):
    """export a csv containing the basin average value

    Args:
        fc (ee.FeatureCollection): e.g. basins
        image (ee.Image): e.g. DEM

    Returns:
        DataFrame: it has fields of 'STAID' and 'mean'
    """
    import ee
    import pandas as pd

    # export as csv
    reducer = kwargs.pop('reducer', ee.Reducer.mean())
    scale = kwargs.pop('scale', image.projection().nominalScale())

    if image.projection().crs().getInfo() != 'EPSG:4326':
        image = image.reproject('EPSG:4326', None, scale)

    means = image.reduceRegions(fc, reducer=reducer, scale=scale, **kwargs)

    if to_drive is None:
        url = means.select(fields, retainGeometry=False).getDownloadURL(filetype='csv')
        if return_url:
            return url
        else:
            return pd.read_csv(url)
    else:
        ee.batch.Export.table.toDrive(
            means, description=to_drive, selectors=fields
        ).start()


def gee_soilgrids(band):
    """download and preprocess soilgrids from google earth engine

    https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/access_on_gee.md

    Args:
        band (string): one of ['bdod', 'cfvo', 'clay', 'sand', 'silt']

    Returns:
        image: ee.Image contains specific band
    """
    import ee

    image = ee.Image('projects/soilgrids-isric/' + band + '_mean')

    # calculate depth-weighted value
    RAW_NAMES = ['B5', 'B10', 'B15', 'B30', 'B40', 'B100']
    RAW_BANDS = [
        '0-5cm_mean',
        '5-15cm_mean',
        '15-30cm_mean',
        '30-60cm_mean',
        '60-100cm_mean',
        '100-200cm_mean',
    ]
    RAW_BANDS = [band + '_' + i for i in RAW_BANDS]
    RAW_DICT = {k: image.select(v) for k, v in zip(RAW_NAMES, RAW_BANDS)}
    image = image.expression(
        '5 * B5 + 10 * B10 + 15 * B15 + 30 * B30 + 40 * B40 + 100 * B100', RAW_DICT
    ).divide(200)
    image = image.reduceResolution(
        reducer=ee.Reducer.mean(), bestEffort=True, maxPixels=64
    )

    return image


def gee_wind(image, name='wind_10m'):
    wind_10m = image.expression(
        'sqrt(u**2 + v**2)',
        {
            'u': image.select('u_component_of_wind_10m'),
            'v': image.select('v_component_of_wind_10m'),
        },
    ).rename(name)
    time = image.get('system:time_start')
    return wind_10m.set('system:time_start', time)


def gee_group_by_month(images):
    # Group by month, and then reduce within groups by mean()
    # the result is an ImageCollection with one image for each
    # month.
    import ee

    months = ee.List.sequence(1, 12)
    by_month = ee.ImageCollection.fromImages(
        months.map(
            lambda m: images.filter(ee.Filter.calendarRange(m, m, 'month'))
            .mean()
            .set('month', m)
        )
    )
    return by_month


def gee_seasonality_index(images):
    avr = images.mean()
    SI = (
        gee_group_by_month(images)
        .map(lambda im: im.subtract(avr).abs())
        .sum()
        .divide(avr.multiply(12))
    )
    return SI
