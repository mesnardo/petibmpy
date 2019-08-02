"""Module to read/write a PetIBM field variable."""

import h5py
import numpy
from scipy import interpolate


def read_field_hdf5(filepath, name):
    """Read a field from HDF5 file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the HDF5 file.
    name : string
        Name of the field variable.

    Returns
    -------
    field : numpy.ndarray
        The PetIBM field variable as a NumPy array of floats.

    """
    f = h5py.File(str(filepath), 'r')
    field = f[name][:]
    f.close()
    return field


def write_field_hdf5(filepath, name, field):
    """Write a field to a HDF5 file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the HDF5 file.
    name : string
        Name of the field variable.
    field : numpy.ndarray
        The PetIBM field variable as a NumPy array of floats.

    """
    f = h5py.File(str(filepath), 'w')
    f.create_dataset(name, data=field)
    f.close()
    return


def linear_interpolation(u, x, xi):
    """Perform a linear interpolation along the first axis.

    Parameters
    ----------
    u : numpy.ndarray
        Array to interpolate.
    x : numpy.ndarray
        Gridline locations.
    xi : float
        Target location.

    Returns
    -------
    ui : numpy.ndarray or float
        Interpolated values.

    """
    idx = numpy.where(x < xi)[0][-1]
    u0, u1 = u[idx], u[idx + 1]
    x0, x1 = x[idx], x[idx + 1]
    xd = (xi - x0) / (x1 - x0)
    ui = (1 - xd) * u0 + xd * u1
    return ui


def interpolate3d(field, grid1, grid2, **kwargs):
    """Interpolate a 3D field from one grid to another.

    Parameters
    ----------
    field : numpy.ndarray
        The 3D field to interpolate.
    grid1 : tuple of numpy.ndarray objects
        The grid on which the field is defined.
        The grid should be provided as (x, y, z).
    grid2 : tuple of numpy.ndarray objects
        The grid on which to interpolate the field.
        The grid should be provided as (x, y, z).
    **kwargs : Arbitrary keyword arguments
        To be passed to scipy.interpolate.interpn.

    Returns
    -------
    field2 : numpy.ndarray
        The interpolated 3D field.

    """
    x1, y1, z1 = grid1
    x2, y2, z2 = grid2
    n2 = x2.size * y2.size * z2.size
    grid = numpy.array(numpy.meshgrid(z2, y2, x2, indexing='ij'))
    grid = numpy.rollaxis(grid, 0, 4).reshape(n2, 3)
    field2 = interpolate.interpn((z1, y1, x1), field, grid, **kwargs)
    field2 = field2.reshape((z2.size, y2.size, x2.size))
    return field2


def interpolate2d(field, grid1, grid2, **kwargs):
    """Interpolate a 2D field from one grid to another.

    Parameters
    ----------
    field : numpy.ndarray
        The 2D field to interpolate.
    grid1 : tuple of numpy.ndarray objects
        The grid on which the field is defined.
        The grid should be provided as (x, y).
    grid2 : tuple of numpy.ndarray objects
        The grid on which to interpolate the field.
        The grid should be provided as (x, y).
    **kwargs : Arbitrary keyword arguments
        To be passed to scipy.interpolate.interpn.

    Returns
    -------
    field2 : numpy.ndarray
        The interpolated 2D field.

    """
    x1, y1 = grid1
    x2, y2 = grid2
    n2 = x2.size * y2.size
    grid = numpy.array(numpy.meshgrid(y2, x2, indexing='ij'))
    grid = numpy.rollaxis(grid, 0, 3).reshape(n2, 2)
    field2 = interpolate.interpn((y1, x1), field, grid, **kwargs)
    field2 = field2.reshape((y2.size, x2.size))
    return field2
