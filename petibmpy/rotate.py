"""Module with function to rotate a geometry."""

import math
import numpy


def rotation2d(x, y, center=(0.0, 0.0), angle=0.0, mode='deg'):
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


def rotation3d(x, y, z,
               roll=0.0, yaw=0.0, pitch=0.0,
               center=(0.0, 0.0, 0.0)):
    """Rotate 3D point.

    Parameters
    ----------
    x : float
        x-coordinate of point.
    y : float
        y-coordinate of point.
    z : float
        z-coordinate of point.
    roll : float (optional)
        Roll angle (in radians); default: 0.0.
    yaw : float (optional)
        Yaw angle (in radians); default: 0.0.
    pitch : float (optional)
        Pitch angle (in radians); default: 0.0.
    center : tuple of floats
        Coordinates of the center of rotation;
        default: [0.0, 0.0, 0.0].

    Returns
    -------
    xr : float
        x-coordinate of rotated point.
    yr : float
        y-coordinate of rotated point.
    zr : float
        z-coordinate of rotated point.

    """
    center = numpy.array(center)
    Rx = numpy.array([[1.0, 0.0, 0.0],
                      [0.0, math.cos(roll), math.sin(roll)],
                      [0.0, -math.sin(roll), math.cos(roll)]])
    Ry = numpy.array([[math.cos(yaw), 0.0, math.sin(yaw)],
                      [0.0, 1.0, 0.0],
                      [-math.sin(yaw), 0.0, math.cos(yaw)]])
    Rz = numpy.array([[math.cos(pitch), math.sin(pitch), 0.0],
                      [-math.sin(pitch), math.cos(pitch), 0.0],
                      [0.0, 0.0, 1.0]])
    point = numpy.array([x, y, z])
    new = Rx.dot(Ry.dot(Rz.dot(point - center))) + center
    xr, yr, zr = new
    return xr, yr, zr


rotation3d_vec = numpy.vectorize(rotation3d,
                                 excluded=['roll', 'yaw', 'pitch',
                                           'center'])
