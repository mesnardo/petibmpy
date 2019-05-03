"""Module to read/write a PetIBM field variable."""

import h5py
import numpy


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
