import numpy as np
from geospace.utils import ds_name


def grid_area(lon_res, lat_res, lat_upper, n_lat):
    """Calculates the area of grid cells on a lat/lon grid in square meters.

    This function implements the equation from Santini et al. 2010 to calculate
    the area of a grid cell on the WGS-84 ellipsoid. The formula is:

    S = (λ_2 - λ_1)(sinφ_2 - sinφ_1)R^2

    where:
    - S is the surface area of the cell.
    - λ_1, λ_2 are the authalic longitudes in radians.
    - φ_1, φ_2 are the authalic latitudes in radians.
    - R is the Earth's authalic radius.

    The authalic radius is the radius of a perfect sphere with the same surface
    area as the reference ellipsoid.

    Reference:
    Santini, M., Taramelli, A., & Sorichetta, A. (2010). ASPHAA: A GIS-Based
    Algorithm to Calculate Cell Area on a Latitude-Longitude (Geographic)
    Regular Grid. Transactions in GIS.
    https://doi.org/10.1111/j.1467-9671.2010.01200.x

    Args:
        lon_res (float): The longitude resolution in degrees.
        lat_res (float): The latitude resolution in degrees.
        lat_upper (float): The upper geodetic latitude of the grid.
        n_lat (int): The number of latitude rows.

    Returns:
        np.ndarray: An array containing the area of each grid cell in square meters.
    """

    # semi-major axis and flattening of the WGS 84
    # https://en.wikipedia.org/wiki/World_Geodetic_System
    SMA = 6378137
    f = 1 / 298.257223563
    e = np.sqrt(2 * f - f * f)
    e1 = 1 - e * e
    e2 = 1 / (2 * e)
    ee = e * e
    sin90 = np.sin(np.pi / 2)
    q = e1 * (
        sin90 / (1 - ee * sin90 * sin90)
        - e2 * np.log((1 - e * sin90) / (1 + e * sin90))
    )
    R2 = SMA * SMA * q / 2

    lat_below = lat_upper - lat_res * n_lat
    sin_lats = np.sin(np.radians(np.linspace(lat_upper, lat_below, n_lat + 1)))
    q_lats = e1 * (
        sin_lats / (1 - ee * sin_lats * sin_lats)
        - e2 * np.log((1 - e * sin_lats) / (1 + e * sin_lats))
    )
    return R2 * np.radians(lon_res) * np.abs((q_lats[:-1] - q_lats[1:]) / q)


def area_per_row(trans, proj_wkt, n_rows, offset=0):
    """Calculates the area of each row in a raster dataset.

    Args:
        trans (tuple): The geotransform of the raster dataset.
        proj_wkt (str): The projection of the raster dataset in WKT format.
        n_rows (int): The number of rows in the raster dataset.
        offset (int, optional): The row offset. Defaults to 0.

    Returns:
        np.ndarray: An array containing the area of each row in square kilometers.
    """
    if 'PROJCS' in proj_wkt:
        cell_area = abs(trans[1] * trans[5]) / 1000000
        return np.full(n_rows, cell_area)

    return (
        grid_area(abs(trans[1]), abs(trans[5]), trans[3] + offset * trans[5], n_rows)
        / 1000000
    )


def real_area(ds, rows, offset=None, return_row_area=False):
    """Calculates the total area of a set of rows in a raster dataset.

    Args:
        ds (gdal.Dataset or str): The raster dataset or its path.
        rows (np.ndarray): An array of row indices.
        offset (int, optional): The row offset. If None, it is calculated from the
                                minimum row index. Defaults to None.
        return_row_area (bool, optional): If True, return the area of each row.
                                          Defaults to False.

    Returns:
        float or np.ndarray: The total area in square kilometers, or an array of row areas.
    """
    rows = np.array(rows).flatten()
    if len(rows) == 0:
        return np.array([])
    ds = ds_name(ds)[0]
    trans = ds.GetGeoTransform()
    proj = ds.GetProjection()
    if 'PROJCS' in proj:
        cell_area = abs(trans[1] * trans[5]) / 1000000
        if return_row_area:
            return np.full(rows.shape, cell_area)
        return rows.shape[0] * cell_area

    if offset is None:
        offset = np.min(rows)
        rows = rows - offset
    row_counts = np.bincount(rows)
    row_area = (
        grid_area(
            abs(trans[1]),
            abs(trans[5]),
            trans[3] + offset * trans[5],
            row_counts.shape[0],
        )
        / 1000000
    )
    if return_row_area:
        return row_area[rows]
    return (row_counts * row_area).sum()


def distance(a, b):
    """Calculates the distance between two points on the WGS-84 ellipsoid.

    This function uses Vincenty's formula to calculate the distance.

    Args:
        a (tuple): The longitude and latitude of the first point (lon, lat).
        b (tuple): The longitude and latitude of the second point (lon, lat).

    Returns:
        float: The distance between the two points in kilometers.
    """
    import math

    ELLIPSOIDS = {
        # model           major (km)   minor (km)     flattening
        'WGS-84': (6378.137, 6356.7523142, 1 / 298.257223563),
    }

    lng1, lat1 = math.radians(a[0]), math.radians(a[1])
    lng2, lat2 = math.radians(b[0]), math.radians(b[1])

    major, minor, f = ELLIPSOIDS['WGS-84']

    delta_lng = lng2 - lng1

    reduced_lat1 = math.atan((1 - f) * math.tan(lat1))
    reduced_lat2 = math.atan((1 - f) * math.tan(lat2))

    sin_reduced1, cos_reduced1 = math.sin(reduced_lat1), math.cos(reduced_lat1)
    sin_reduced2, cos_reduced2 = math.sin(reduced_lat2), math.cos(reduced_lat2)

    lambda_lng = delta_lng
    lambda_prime = 2 * math.pi

    iter_limit = 20

    i = 0

    while i == 0 or (abs(lambda_lng - lambda_prime) > 10e-12 and i <= iter_limit):
        i += 1

        sin_lambda_lng, cos_lambda_lng = math.sin(lambda_lng), math.cos(lambda_lng)

        sin_sigma = math.sqrt(
            (cos_reduced2 * sin_lambda_lng) ** 2
            + (
                cos_reduced1 * sin_reduced2
                - sin_reduced1 * cos_reduced2 * cos_lambda_lng
            )
            ** 2
        )

        if sin_sigma == 0:
            return 0  # Coincident points

        cos_sigma = (
            sin_reduced1 * sin_reduced2 + cos_reduced1 * cos_reduced2 * cos_lambda_lng
        )

        sigma = math.atan2(sin_sigma, cos_sigma)

        sin_alpha = cos_reduced1 * cos_reduced2 * sin_lambda_lng / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2

        if cos_sq_alpha != 0:
            cos2_sigma_m = cos_sigma - 2 * (sin_reduced1 * sin_reduced2 / cos_sq_alpha)
        else:
            cos2_sigma_m = 0.0  # Equatorial line

        C = f / 16.0 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))

        lambda_prime = lambda_lng
        lambda_lng = delta_lng + (1 - C) * f * sin_alpha * (
            sigma
            + C
            * sin_sigma
            * (cos2_sigma_m + C * cos_sigma * (-1 + 2 * cos2_sigma_m**2))
        )

    if i > iter_limit:
        raise ValueError('Vincenty formula failed to converge!')

    u_sq = cos_sq_alpha * (major**2 - minor**2) / minor**2

    A = 1 + u_sq / 16384.0 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))

    B = u_sq / 1024.0 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

    delta_sigma = (
        B
        * sin_sigma
        * (
            cos2_sigma_m
            + B
            / 4.0
            * (
                cos_sigma * (-1 + 2 * cos2_sigma_m**2)
                - B
                / 6.0
                * cos2_sigma_m
                * (-3 + 4 * sin_sigma**2)
                * (-3 + 4 * cos2_sigma_m**2)
            )
        )
    )

    s = minor * A * (sigma - delta_sigma)
    return s


if __name__ == '__main__':
    area_columns = grid_area(1 / 3600, 1 / 3600, 60, 300000)
    area_columns = grid_area(0.5, 0.5, 60, 100)
