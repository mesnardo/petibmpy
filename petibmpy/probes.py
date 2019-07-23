"""Module for the PetIBM probes."""

import h5py
import numpy


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


def get_probe_volume_limits(x, loc, neighbors=1, btype='both', decimals=6):
    """Return the interval limits (left, right, or both).

    Parameters
    ----------
    x : 1D numpy.ndarray
        The gridline coordinates.
    loc : float
        Location for which to request the gridline boundaries.
    neighbors : integer (optional)
        Number of neighbors to include as a buffer; default: 1.
    btype : string (optional)
        Type of boundary to return;
        choices: 'left', 'right, or 'both';
        default: 'both'.
    decimals : integer (optional)
        Number of digits used to round numbers; default: 6.

    Returns
    -------
    'left', 'right', or 'both' gridline boundaries.

    """
    # Check the value of btype.
    if btype not in ['left', 'right', 'both']:
        raise ValueError("Wrong type! "
                         "(Use btype = 'left', 'right', or 'both'.)")
    x = numpy.round(x, decimals=decimals)
    # Get the index of the gridline point just before the point of interest.
    if btype == 'left':
        idx = numpy.where(x > loc)[0][0] - 1
    # Get the index of the gridline point just after the point of interest.
    elif btype == 'right':
        idx = numpy.where(x < loc)[0][-1] + 1
    # Get the index of the gridline point closest to the point of interest.
    else:
        idx = numpy.argmin(numpy.abs(x - loc))
    # Store the left and right boundary points.
    left = 0.5 * (x[idx - neighbors - 1] + x[idx - neighbors])
    right = 0.5 * (x[idx + neighbors] + x[idx + neighbors + 1])
    # Round the values.
    left = round(float(left), ndigits=decimals)
    right = round(float(right), ndigits=decimals)
    # Return the left, right, or both boundary points.
    if btype == 'left':
        return left
    elif btype == 'right':
        return right
    else:
        return (left, right)
