import numpy as np
from scipy.optimize import fsolve, root
from scipy.signal import argrelextrema

def find_inflection_points(x, y):
    """
    Determines the inflection points of the function based on second derivative.

    Parameters:
    - x, y: data points

    Returns:
    - x and y values of turning points and their indices.
    """
    dy = np.diff(y) / np.diff(x)
    dx = (x[:-1] + x[1:]) / 2
    d2y = np.diff(dy) / np.diff(dx)
    sign_changes = np.diff(np.sign(d2y))
    i = np.where(sign_changes)[0] + 2
    return i, x[i], y[i]

def find_common_tangent(f, g, xL_init, xR_init):
    """
    Determines the common tangent for two functions.

    Parameters:
    - f, g: Polynomial functions
    - xL_init, xR_init: Initial guesses for intersection points

    Returns:
    - Intersection points and y-values.
    """
    def fobj(x):
        xL, xR = x
        df = f.deriv()(xL)
        dg = g.deriv()(xR)
        dy = (g(xR) - f(xL)) / (xR - xL)
        return [df - dg, dy - df]
    xL, xR = fsolve(fobj, x0=[xL_init, xR_init])
    # [xL, xR], info, *_ = fsolve(fobj, x0=[xL_init, xR_init], full_output=True)
    # print(info['nfev'])
    # kwargs = dict(method='krylov', options={
    #               'maxiter': 20, 'xtol': 1e-3})
    # res = root(fobj, x0=[xL_init, xR_init], **kwargs)
    # [xL, xR] = res.x
    # print(res.nit)
    return xL, xR, f(xL), g(xR)

def approximate_between_points(x, y, start, end, deg=5):
    segment_x = x[start:end+1]
    segment_y = y[start:end+1]
    midpoint = (segment_x[0] + segment_x[-1]) / 2
    params = np.polyfit(segment_x, segment_y, deg=deg)
    func = np.poly1d(params)
    return func, segment_x, midpoint

def get_segment_border_indices(x, start, end, fraction=0.5, min_values=5):
    # Calculate the intermediate point using the given fraction
    intermediate_x = x[start] + fraction * (x[end] - x[start])
    # Find the index in x that is closest to this intermediate point
    # closest_index = start + np.argmin(np.abs(x[start:end+1] - intermediate_x))
    closest_index = np.argmin(np.abs(x - intermediate_x))
    # Ensure that there are at least min_values between start and closest_index
    if abs(closest_index - start) < min_values:
        closest_index = end
    return sorted([start, closest_index])

def find_binodal(x1, gmix, iL, iR, ax=None):
    # Approximate segments between pure component and inflection points
    borderL = get_segment_border_indices(x1, start=iL, end=0)
    borderR = get_segment_border_indices(x1, start=iR, end=-1)
    f, xL_range, mL = approximate_between_points(x1, gmix, *borderL, deg=5)
    g, xR_range, mR = approximate_between_points(x1, gmix, *borderR, deg=5)
    # Find common tangent
    xL, xR, yL, yR = find_common_tangent(f, g, mL, mR)

    # Check if results are outside the fitting range
    adjustL = xL < xL_range[0] or xL > xL_range[-1]
    adjustR = xR < xR_range[0] or xR > xR_range[-1]

    # If outside, adjust the respective range and recalculate xL and xR
    if adjustL or adjustR:
        if adjustL:
            borderL = get_segment_border_indices(
                x1, start=iL, end=0, fraction=1)
            f, xL_range, mL = approximate_between_points(
                x1, gmix, *borderL, deg=5)
        if adjustR:
            borderR = get_segment_border_indices(
                x1, start=iR, end=-1, fraction=1)
            g, xR_range, mR = approximate_between_points(
                x1, gmix, *borderR, deg=10)

        # Find common tangent
        xL, xR, yL, yR = find_common_tangent(f, g, mL, mR)

    if ax is not None:
        ax.plot(xL_range, f(xL_range), 'r', lw=1)
        ax.plot(xR_range, g(xR_range), 'r', lw=1)

    return xL, xR, yL, yR

def find_local_extremum(y, typ='minimum'):
    kind = np.less if "min" in typ else np.greater
    return argrelextrema(y, kind)[0]

def approx_binodal(x1, gmix, xS, iL, iR):
    x, y = x1, gmix
    local_minima = find_local_extremum(y, 'minimum')
    min_x = x[local_minima]

    # Initialize return values
    xL = xR = yL = yR = None

    if len(min_x) == 0:
        pass
    # Check the number of local minima found
    if len(min_x) == 1:
        # Precompute the boolean condition to avoid redundancy
        is_greater_than_mid_x = min_x > max(xS)
        # Select indices based on the precomputed condition
        i = iL if is_greater_than_mid_x else iR
        slice_ = slice(None, i+1) if is_greater_than_mid_x else slice(i, None)

        # Slicing the arrays based on the selected indices
        x1, y1 = x[slice_], y[slice_]

        # Fit a line to the boundary points
        m, n = np.polyfit(x1[[0, -1]], y1[[0, -1]], 1)
        f = np.poly1d([m, n])

        # Calculate the difference for finding the local minima
        diff = y1 - f(x1)
        local_minima = find_local_extremum(diff, 'minimum')

        # Adjust local_minima index based on the original x array
        if not is_greater_than_mid_x:
            local_minima += i  # Adjust index for the slice offset

        # Check if new local minima were found and sort
        if local_minima.size > 0:
            # Sort to ensure xL < xR
            xL, xR = np.sort([min_x[0], x[local_minima][0]])
            yL, yR = y[np.where((x == xL) | (x == xR))]

    elif len(min_x) == 2:
        xL, xR = min_x
        yL, yR = y[local_minima]

    else:
        pass

    return xL, xR, yL, yR


def estimate_lle_from_gmix(x1, gmix, rough=True, ax=None):
    # Initialize return values
    xL = xR = yL = yR = None

    # Find inflection points
    try:
        idx, xS, yS = find_inflection_points(x1, gmix)
        iL, iR = idx[:2]
    except ValueError:
        # print("Warning: No inflection points found for current iteration.")
        return xL, xR, yL, yR

    # Find binodal
    if rough:
        xL, xR, yL, yR = approx_binodal(x1, gmix, xS, iL, iR)
    else:
        xL, xR, yL, yR = find_binodal(x1, gmix, iL, iR, ax=ax)

    # Round xL and xR if they are not None
    xL = round(xL, 8) if xL is not None else None
    xR = round(xR, 8) if xR is not None else None

    if ax is not None and None not in (xL, xR, yL, yR):
        # Plot common tangent
        m = (yR-yL)/(xR-xL)
        n = yR - m*xR
        t = np.poly1d([m, n])
        ylim = ax.get_ylim()
        ax.plot([xL, xR], [yL, yR], 'ko', mfc='w', zorder=3)
        ax.plot(x1, t(x1), 'k:')
        ax.set_ylim(ylim)

    return xL, xR, yL, yR
