"""Module for the PetIBM probes."""

import h5py


def read_probe_volume_hdf5(filepath, name, time):
    """Read a volume field from file at a given recorded time.

    Parameters
    ----------
    filepath : pathlib.Path or string
        Path of the file to read from.
    name : string
        Name of the field to read.
    time : float
        Time value at which the volume field was recorded.

    Returns
    -------
    mesh : tuple of numpy.ndarray objects
        The gridline coordinates of the volume.
    values : numpy.ndarray
        The volume field.

    """
    f = h5py.File(str(filepath), 'r')
    mesh_group = f['mesh']
    dim = len(mesh_group)
    mesh = mesh_group['x'][:], mesh_group['y'][:]
    if dim == 3:
        mesh = mesh + (mesh_group['z'][:],)
    mesh_sizes = [line.size for line in mesh]
    field_group = f[name]
    time_str = '{:0.6f}'.format(round(time, ndigits=6))
    values = field_group[time_str][:]
    values = values.reshape(mesh_sizes[::-1])
    f.close()
    return mesh, values
