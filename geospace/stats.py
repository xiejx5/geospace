import numpy as np


def mean(arr, weights):
    """
    Calculates the area-weighted mean.

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).
        weights (np.ndarray): A 2D array of weights with shape (lat, lon).

    Returns:
        np.ndarray: The weighted mean for each band.
    """

    return np.ma.average(arr, weights=weights, axis=(1, 2)).filled(np.nan)


def max(arr):
    """
    Calculates the maximum value (unweighted).

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).

    Returns:
        np.ndarray: The maximum value for each band.
    """

    result = np.ma.max(arr, axis=(1, 2))
    if np.issubdtype(result.dtype, np.integer):
        result = result.astype(float)
    return result.filled(np.nan)


def min(arr):
    """
    Calculates the minimum value (unweighted).

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).

    Returns:
        np.ndarray: The minimum value for each band.
    """

    result = np.ma.min(arr, axis=(1, 2))
    if np.issubdtype(result.dtype, np.integer):
        result = result.astype(float)
    return result.filled(np.nan)


def median(arr):
    """
    Calculates the median value (unweighted).

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).

    Returns:
        np.ndarray: The median value for each band.
    """

    result = np.ma.median(arr, axis=(1, 2))
    if np.issubdtype(result.dtype, np.integer):
        result = result.astype(float)
    return result.filled(np.nan)


def sum(arr, weights):
    """
    Calculates the area-weighted sum (total value * area).

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).
        weights (np.ndarray): A 2D array of weights with shape (lat, lon).

    Returns:
        np.ndarray: The weighted sum for each band.
    """

    return np.ma.sum(arr * weights, axis=(1, 2), dtype=float).filled(np.nan)


def area(arr, weights):
    """
    Calculates the total area of intersected valid pixels.

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).
        weights (np.ndarray): A 2D array of weights with shape (lat, lon).

    Returns:
        np.ndarray: The sum of weights for valid pixels for each band.
    """

    # Create a mask for weights based on arr mask
    masked_weights = np.ma.array(np.broadcast_to(weights, arr.shape), mask=arr.mask)
    return np.ma.sum(masked_weights, axis=(1, 2), dtype=float).filled(np.nan)


def std(arr, weights):
    """
    Calculates the area-weighted standard deviation.

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).
        weights (np.ndarray): A 2D array of weights with shape (lat, lon).

    Returns:
        np.ndarray: The weighted standard deviation for each band.
    """

    avg = np.ma.average(arr, weights=weights, axis=(1, 2))
    variance = np.ma.average(
        (arr - avg[:, None, None]) ** 2, weights=weights, axis=(1, 2)
    )
    return np.sqrt(variance).filled(np.nan)


def mode(arr):
    """
    Calculates the mode (most frequent value) (unweighted).

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).

    Returns:
        np.ndarray: The mode for each band.
    """
    arr = arr.reshape(arr.shape[0], -1)

    def _get_mode(a):
        # working with masked array 1d
        compressed = a.compressed()
        if compressed.size == 0:
            return np.nan
        values, counts = np.unique(compressed, return_counts=True)
        return values[np.argmax(counts)]

    return np.apply_along_axis(_get_mode, 1, arr)


def _get_max_indices(arr, weights):
    """
    Helper to find the indices (y, x) of the pixel with the largest weight for each band.

    Returns:
        tuple: (y, x)
            y (np.ndarray): Row indices of max weight pixels (float, with NaNs for masked).
            x (np.ndarray): Column indices of max weight pixels (float, with NaNs for masked).
    """
    # Create a mask (n_band, n_pixel) for weights based on arr mask
    masked_weights = np.ma.array(np.broadcast_to(weights, arr.shape), mask=arr.mask)
    masked_weights = masked_weights.reshape(arr.shape[0], -1)

    # Convert linear max indices to y, x
    max_indices = np.ma.argmax(masked_weights, axis=1)
    y, x = np.unravel_index(max_indices, weights.shape)

    # Convert to float to support NaNs
    y = y.astype(float)
    x = x.astype(float)

    # Set indices to NaN where the band is fully masked
    all_masked = masked_weights.mask.all(axis=1)
    y[all_masked] = np.nan
    x[all_masked] = np.nan

    return y, x


def max_area_lon(arr, weights, trans):
    """
    Calculates the longitude of the pixel centroid with the largest intersected area.

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).
        weights (np.ndarray): A 2D array of weights with shape (lat, lon).
        trans (tuple): The GDAL geotransform tuple (ulx, xres, xskew, uly, yskew, yres).

    Returns:
        np.ndarray: The longitude of the pixel with the largest intersection for each band.
    """
    y, x = _get_max_indices(arr, weights)

    # Calculate coordinates
    lons = trans[0] + (x + 0.5) * trans[1] + (y + 0.5) * trans[2]

    return lons


def max_area_lat(arr, weights, trans):
    """
    Calculates the latitude of the pixel centroid with the largest intersected area.

    Args:
        arr (np.ma.MaskedArray): A 3D masked array with shape (n_bands, lat, lon).
        weights (np.ndarray): A 2D array of weights with shape (lat, lon).
        trans (tuple): The GDAL geotransform tuple (ulx, xres, xskew, uly, yskew, yres).

    Returns:
        np.ndarray: The latitude of the pixel with the largest intersection for each band.
    """
    y, x = _get_max_indices(arr, weights)

    # Calculate coordinates
    lats = trans[3] + (x + 0.5) * trans[4] + (y + 0.5) * trans[5]

    return lats
