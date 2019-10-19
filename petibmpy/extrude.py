"""Module with function to extrude a 2D geometry in the third direction."""

import math
import numpy


def extrude2d(x, y, limits, n=None, ds=None, force=False):
    """Extrude the two-dimensional section along the third direction (z).

    Parameters
    ----------
    x : numpy.ndarray
        x-coordinates of the section.
    y : numpy.ndarray
        y-coordinates of the section.
    limits : 2-list of floats
        Limits of the extrusion.
    n : integer, optional
        Number of divisions in the z-direction;
        default: None.
    ds : float, optional
        Desired segment-length;
        default: None.
    force : boolean, optional
        Forces the extrusion to the limits prescribed;
        default: False.

    Returns
    -------
    x : numpy.ndarray
        x-coordinates of the geometry.
    y : numpy.ndarray
        y-coordinates of the geometry.
    z : numpy.ndarray
        z-coordinates of the geometry.

    """
    if not (ds or n):
        raise ValueError('both ds and n are set to None')
    elif abs(limits[0] - limits[1]) < 1e-6:
        raise ValueError('limits are too close from each other')
    z_start, z_end = limits
    if not n:
        n = int(math.ceil(abs(z_start - z_end) / ds))
    ds = abs(z_start - z_end) / n
    s = math.copysign(1.0, z_end - z_start)
    if force:
        z = numpy.linspace(z_start, z_end, n + 1)
    else:
        z = numpy.linspace(z_start + s * 0.5 * ds, z_end - s * 0.5 * ds, n)
    nx, ny, nz = x.size, y.size, z.size
    x, y = numpy.tile(x, nz), numpy.tile(y, nz)
    z = numpy.tile(z, (nx, 1)).T.flatten()
    return x, y, z
