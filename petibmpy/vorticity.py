"""Module with functions to compute the vorticity."""

import numpy


def gradient(u, grid, axis=0):
    """Compute the gradient of u along a given axis.

    Parameters
    ----------
    u : numpy.ndarray
        An N-dimensional array.
    grid : tuple of N 1-D arrays
        Grid on which u is defined.
    axis : int, optional
        Axis along which to compute the gradient, by default 0.

    Returns
    -------
    numpy.ndarray
        The N-dimensional gradient of u along a given axis.
    tuple of N 1-D arrays
        The grid on which is defined the gradient.

    """
    u = numpy.moveaxis(u, axis, 0)
    grid_N = numpy.meshgrid(*grid)
    x = grid_N[-1 - axis]
    x = numpy.moveaxis(x, axis, 0)
    grad = _gradient(u, x)
    grad = numpy.moveaxis(grad, 0, axis)
    x = grid[-1 - axis]
    x = 0.5 * (x[:-1] + x[1:])
    return grad, x


def _gradient(u, x):
    grad = (u[1:] - u[:-1]) / (x[1:] - x[:-1])
    return grad


def compute_wx(v, w, grid_v, grid_w):
    """Compute the x-component of the vorticity field.

    Parameters
    ----------
    v : numpy.ndarray
        y-component of the velocity field as a N-dimensional array.
    w : numpy.ndarray
        z-component of the velocity field as a N-dimensional array.
    grid_v : tuple of N 1-D arrays
        Grid on which the y-velocity is defined.
    grid_w : tuple of N 1-D arrays
        Grid on which the z-velocity is defined.

    Returns
    -------
    numpy.ndarray
        x-component of the vorticity field.
    tuple of N 1-D arrays
        Grid on which the x-vorticity is defined.

    """
    dim = len(v.shape)
    assert dim == 3, 'Velocity field should be 3D'
    dwdy, y = gradient(w, grid_w, axis=1)
    dvdz, z = gradient(v, grid_v, axis=0)
    wx = dwdy - dvdz
    x = grid_v[0]
    return wx, (x, y, z)


def compute_wy(u, w, grid_u, grid_w):
    """Compute the y-component of the vorticity field.

    Parameters
    ----------
    u : numpy.ndarray
        x-component of the velocity field as a N-dimensional array.
    w : numpy.ndarray
        z-component of the velocity field as a N-dimensional array.
    grid_x : tuple of N 1-D arrays
        Grid on which the x-velocity is defined.
    grid_w : tuple of N 1-D arrays
        Grid on which the z-velocity is defined.

    Returns
    -------
    numpy.ndarray
        y-component of the vorticity field.
    tuple of N 1-D arrays
        Grid on which the y-vorticity is defined.

    """
    dim = len(u.shape)
    assert dim == 3, 'Velocity field should be 3D'
    dudz, z = gradient(u, grid_u, axis=0)
    dwdx, x = gradient(w, grid_w, axis=2)
    wy = dudz - dwdx
    y = grid_u[1]
    return wy, (x, y, z)


def compute_wz(u, v, grid_u, grid_v):
    """Compute the z-component of the vorticity field.

    Parameters
    ----------
    u : numpy.ndarray
        x-component of the velocity field as a N-dimensional array.
    v : numpy.ndarray
        y-component of the velocity field as a N-dimensional array.
    grid_u : tuple of N 1-D arrays
        Grid on which the x-velocity is defined.
    grid_v : tuple of N 1-D arrays
        Grid on which the y-velocity is defined.

    Returns
    -------
    numpy.ndarray
        z-component of the vorticity field.
    tuple of N 1-D arrays
        Grid on which the z-vorticity is defined.

    """
    dim = len(u.shape)
    assert dim == 2 or dim == 3, 'Velocity field should be 2D or 3D'
    dudy, y = gradient(u, grid_u, axis=0 + (dim == 3))
    dvdx, x = gradient(v, grid_v, axis=1 + (dim == 3))
    wz = dvdx - dudy
    if dim == 2:
        return wz, (x, y)
    z = grid_u[-1]
    return wz, (x, y, z)
