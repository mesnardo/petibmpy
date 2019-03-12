"""Module with function to rotate a 2D geometry."""

import numpy


def rotate2d(x, y, center=(0.0, 0.0), angle=0.0, mode='deg'):
    """Rotate (x, y) coordinates around a center.

    Parameters
    ----------
    x : numpy.ndarray of floats
        The x-coordinates to rotate.
    y : numpy.ndarray of floats
        The y-coordinates to rotate.
    center : 2-tuple of floats, optional
        The center of rotation; default: (0.0, 0.0).
    angle : float, optional
        The angle of rotation; default: 0.0.
    mode : string, optional
        Whether angle is provided in degrees or in radians;
        choices: 'deg', 'rad'; default: 'deg'.

    Returns
    -------
    x_new : numpy.ndarray of floats
        The rotated x-coordinates.
    y_new : numpy.ndarray of floats
        The rotated y-coordinates.

    """
    if mode not in ['deg', 'rad']:
        raise ValueError('mode should set to either "deg" or "rad"')
    if mode == 'deg':
        angle = numpy.radians(angle)
    xc, yc = center
    x_new = xc + (x - xc) * numpy.cos(angle) - (y - yc) * numpy.sin(angle)
    y_new = yc + (x - xc) * numpy.sin(angle) + (y - yc) * numpy.cos(angle)
    return x_new, y_new
