"""Module to read/write a PetIBM grid."""

import h5py


def read_grid_hdf5(filepath, name):
    """Read a grid from HDF5 file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the HDF5 file.
    name : string
        Name of the grid.

    Returns
    -------
    x : numpy.ndarray
        The x-coordinates along a gridline in the x-direction.
    y : numpy.ndarray
        The y-coordinates along a gridline in the y-direction.
    z : numpy.ndarray
        The z-coordinates along a gridline in the z-direction.

    """
    f = h5py.File(str(filepath), 'r')
    dim = len(f[name])
    x, y, z = f[name]['x'][:], f[name]['y'][:], None
    if dim == 3:
        z = f[name]['z'][:]
    f.close()
    if z is None or len(z) == 1:
        return x, y
    return x, y, z


def write_grid_hdf5(filepath, name, *grid):
    """Write a grid to a HDF5 file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the HDF5 file.
    name : string
        Name of the grid.
    grid : tuple of numpy.ndarray objects
        The gridline coordinates as 1D arrays of floats.

    """
    labels = ('x', 'y', 'z')
    f = h5py.File(str(filepath), 'w')
    group = f.create_group(name)
    for i, gridline in enumerate(grid):
        group.create_dataset(labels[i], data=gridline)
    f.close()
    return
