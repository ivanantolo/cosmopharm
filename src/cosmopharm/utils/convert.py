import numpy as np


def convert(x: np.ndarray, Mw: np.ndarray, to: str = 'weight') -> np.ndarray:
    """
    Convert between mole fraction (x) and mass/weight fraction (w) for a mixture.

    Parameters:
        x (numpy.ndarray): A 1D or 2D NumPy array representing mole fractions of components in the mixture.
                           - If 1D, it should have the shape (n,) with 'n' as the number of components.
                           - If 2D, it should have the shape (n, m) where 'n' is the number of components
                             and 'm' is the number of data points.
        to (str, optional): The conversion direction. Should be either 'weight' (default) to convert from mole
                           fraction to weight fraction, or 'mole' to convert from weight fraction to mole fraction.

    Returns:
        numpy.ndarray: The converted fractions, with the same shape as 'x'.

    Example:
        >>> sle = SLE(mix=mix)  # Replace with the actual mixture object
        >>> x = np.array([0.4, 0.6])  # Mole fractions of two components
        >>> w = mix.convert(x, to='weight')  # Convert to weight fractions
        >>> print(w)
        array([[0.01373165],
               [0.98626835]])
    """

    # Check if x is a NumPy array
    if not isinstance(x, np.ndarray):
        raise ValueError("Input 'x' must be a 1D or 2D NumPy array.")

    # Check if x is a scalar (0-dimensional)
    if x.shape == ():
        raise ValueError(
            "Input 'x' must be a 1D or 2D NumPy array, not a scalar.")

    if to not in ['weight', 'mole']:
        raise ValueError("Invalid 'to' argument. Use 'weight' or 'mole'.")

    # Check and reshape Mw if needed
    Mw = Mw[:, np.newaxis] if Mw.ndim == 1 else Mw  # Reshape to (n, 1)
    # Check and reshape x if needed
    x = x[:, np.newaxis] if x.ndim == 1 else x  # Reshape to (n, 1)

    # Check if Mw and x have the same number of components
    if Mw.shape[0] != x.shape[0]:
        raise ValueError(
            "Number of components in 'Mw' and 'x' must match.")

    # Calculate the numerator for conversion based on 'to' argument
    num = x * Mw if to == 'weight' else x / Mw
    # Calculate the denominator for conversion
    den = np.sum(num, axis=0)
    # Calculate the final conversion using the numerator and denominator
    res = num / den
    # Replace nan values with 0
    res = np.nan_to_num(res)
    return res
