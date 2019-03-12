"""Module to read/write a PetIBM field variable."""

import h5py


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
