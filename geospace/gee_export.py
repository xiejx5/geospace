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
    """Initializes the Google Earth Engine API.

    This function checks for an internet connection and sets a proxy if needed.
    It then authenticates and initializes the Earth Engine API.

    Args:
        proxy (str, optional): The proxy URL to use. Defaults to 'http://127.0.0.1:7890'.
    """
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
    """Exports a Google Earth Engine image as a GeoTIFF file.

    This function downloads the image directly by generating a download URL.
    It is suitable for smaller images. For larger images, consider using
    `gee_to_drive`.

    Args:
        image (ee.Image): The Earth Engine image to export.
        filename (str): The path to the output GeoTIFF file.
        timeout (int, optional): The timeout for the download in seconds. Defaults to 300.
        **params: Additional parameters for ee.Image.getDownloadURL().
            crs (str, optional): A CRS string to use for bands that do not specify one.
            crs_transform (list, optional): A list of 6 numbers to control the spatial
                resolution and alignment, e.g., [1, 0, -180, 0, -1, 90].
            region (ee.Geometry, list, or str, optional): A polygon specifying a region
                to download. Can be set to 'global' to use the image's full geometry.
            filePerBand (bool, optional): If True, produces a different GeoTIFF per band.
                Defaults to False.
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
    """Exports a Google Earth Engine image to Google Drive.

    This is a wrapper for ee.batch.Export.image.toDrive() that provides
    shortcuts for common export configurations.

    Args:
        image (ee.Image): The image to be exported.
        **params: Keyword arguments for ee.batch.Export.image.toDrive().
            description (str, optional): Human-readable name of the task.
                Defaults to the first band name of the image.
            folder (str, optional): The name of a folder in your Drive account.
            dimensions (str or int, optional): The dimensions of the exported image.
                Can be a single integer for max dimension or "WIDTHxHEIGHT".
                If region is 'global' or 'nonATA', can be set to 'default' for
                automatic dimensions.
            region (list, str, or ee.Geometry, optional): The region to export.
                Can be specified as coordinates. Special string values are also
                accepted:
                - 'global': Exports the entire globe ([-180, -90, 180, 90]).
                - 'nonATA': Exports the globe excluding Antarctica ([-180, -60, 180, 89]).
                Defaults to the image's region.
            crs (str, optional): The coordinate reference system. Defaults to 'EPSG:4326'.
            crsTransform (list or str, optional): The affine transform of the CRS.
                Can be specified as a list of 6 numbers. A special string value
                is also accepted:
                - 'ERA5_LAND': Uses the transform [0.1, 0, -180.05, 0, -0.1, 90.05].
                Defaults to the image's native transform.
            maxPixels (int, optional): The maximum number of pixels to export.
                Defaults to 1e13.
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
    fc, image, columns=['STAID', '.*mean'], return_url=False, to_drive=None, **kwargs
):
    """Exports zonal statistics from a Google Earth Engine image to a CSV file.

    Args:
        fc (ee.FeatureCollection): The feature collection to calculate statistics for (e.g., basins).
        image (ee.Image): The image to extract values from (e.g., DEM).
        columns (list, optional): The columns to include in the output.
                                  Defaults to ['STAID', '.*mean'].
        return_url (bool, optional): If True, return the download URL instead of a DataFrame.
                                     Defaults to False.
        to_drive (str, optional): If not None, export the table to Google Drive with this
                                  description. Defaults to None.
        **kwargs: Additional parameters for ee.Image.reduceRegions().

    Returns:
        pandas.DataFrame or str: A DataFrame with the zonal statistics or the download URL.
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
        url = means.select(columns, retainGeometry=False).getDownloadURL(filetype='csv')
        if return_url:
            return url
        else:
            return pd.read_csv(url)
    else:
        ee.batch.Export.table.toDrive(
            means, description=to_drive, selectors=columns
        ).start()


def gee_soilgrids(band):
    """Downloads and preprocesses SoilGrids data from Google Earth Engine.

    This function calculates a depth-weighted average of the specified soil property.
    https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/access_on_gee.md

    Args:
        band (str): The SoilGrids band to download. One of ['bdod', 'cfvo', 'clay', 'sand', 'silt'].

    Returns:
        ee.Image: An Earth Engine image containing the preprocessed SoilGrids data.
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
    """Calculates wind speed from u and v components of wind.

    Args:
        image (ee.Image): An Earth Engine image with 'u_component_of_wind_10m' and
                          'v_component_of_wind_10m' bands.
        name (str, optional): The name of the output wind speed band. Defaults to 'wind_10m'.

    Returns:
        ee.Image: An Earth Engine image with the calculated wind speed.
    """
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
    """Groups an ImageCollection by month and calculates the mean for each month.

    Args:
        images (ee.ImageCollection): The input ImageCollection.

    Returns:
        ee.ImageCollection: An ImageCollection with one image for each month.
    """
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
    """Calculates the seasonality index of an ImageCollection.

    Args:
        images (ee.ImageCollection): The input ImageCollection.

    Returns:
        ee.Image: An image representing the seasonality index.
    """
    avr = images.mean()
    SI = (
        gee_group_by_month(images)
        .map(lambda im: im.subtract(avr).abs())
        .sum()
        .divide(avr.multiply(12))
    )
    return SI
