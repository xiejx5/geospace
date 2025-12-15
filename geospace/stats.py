import numpy as np


def mean(arr, weights):
    """
    Calculates the area-weighted mean.

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): An array of weights (e.g., area fractions).
                              Shape should match the flattened spatial dimensions of arr.

    Returns:
        np.ndarray: The weighted mean for each band.
    """
    return np.ma.average(arr, weights=weights, axis=1).filled(np.nan)


def max(arr, weights):
    """
    Calculates the maximum value (unweighted).

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): Ignored.

    Returns:
        np.ndarray: The maximum value for each band.
    """
    result = np.ma.max(arr, axis=1)
    if np.issubdtype(result.dtype, np.integer):
        result = result.astype(float)
    return result.filled(np.nan)


def min(arr, weights):
    """
    Calculates the minimum value (unweighted).

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): Ignored.

    Returns:
        np.ndarray: The minimum value for each band.
    """
    result = np.ma.min(arr, axis=1)
    if np.issubdtype(result.dtype, np.integer):
        result = result.astype(float)
    return result.filled(np.nan)


def median(arr, weights):
    """
    Calculates the median value (unweighted).

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): Ignored.

    Returns:
        np.ndarray: The median value for each band.
    """
    result = np.ma.median(arr, axis=1)
    # Median is usually float anyway if calculated, but numpy might verify
    if np.issubdtype(result.dtype, np.integer):
        result = result.astype(float)
    return result.filled(np.nan)


def sum(arr, weights):
    """
    Calculates the area-weighted sum (total value * area).

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): An array of weights (e.g., area fractions).

    Returns:
        np.ndarray: The weighted sum for each band.
    """
    return np.ma.sum(arr * weights, axis=1, dtype=float).filled(np.nan)


def count(arr, weights):
    """
    Calculates the total weighted count (total area).

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): An array of weights (e.g., area fractions).

    Returns:
        np.ndarray: The sum of weights for valid pixels for each band.
    """
    # Create a mask for weights based on arr mask
    masked_weights = np.ma.array(np.broadcast_to(weights, arr.shape), mask=arr.mask)
    return np.ma.sum(masked_weights, axis=1, dtype=float).filled(np.nan)


def std(arr, weights):
    """
    Calculates the area-weighted standard deviation.

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): An array of weights (e.g., area fractions).

    Returns:
        np.ndarray: The weighted standard deviation for each band.
    """
    avg = np.ma.average(arr, weights=weights, axis=1)
    variance = np.ma.average((arr - avg[:, None]) ** 2, weights=weights, axis=1)
    return np.sqrt(variance).filled(np.nan)


def mode(arr, weights):
    """
    Calculates the mode (most frequent value) (unweighted).

    Args:
        arr (np.ma.MaskedArray): A masked array of values.
                                 Shape should be (n_bands, n_pixels).
        weights (np.ndarray): Ignored.

    Returns:
        np.ndarray: The mode for each band.
    """

    def _get_mode(a):
        # working with masked array 1d
        compressed = a.compressed()
        if compressed.size == 0:
            return np.nan
        values, counts = np.unique(compressed, return_counts=True)
        return values[np.argmax(counts)]

    return np.apply_along_axis(_get_mode, 1, arr)
