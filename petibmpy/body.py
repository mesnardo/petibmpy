"""Module with I/O functions for immersed body."""

import numpy


def write_body(filepath, *coords):
    """Save the boundary coordinates to a file.

    Parameters
    ----------
    filepath : pathlib.Path object or string
        Path of the file to write.
    coords : tuple of lists or numpy.ndarray objects
        The x, y, and z coordinates of the boundary.

    """
    with open(filepath, 'w') as outfile:
        outfile.write('{}\n'.format(len(coords[0])))
    with open(filepath, 'ab') as outfile:
        numpy.savetxt(outfile, numpy.c_[coords])


def read_body(filepath, **kwargs):
    """Read the boundary coordinates from a file.

    Parameters
    ----------
    filepath : pathlib.Path object or string
        Path of the file to read.
    kwargs : dictionary
        Keyword arguments to pass to numpy.loadtxt.

    Returns
    -------
    coords : numpy.ndarray
        The boundary coordinates.

    """
    with open(filepath, 'r') as infile:
        coords = numpy.loadtxt(infile, unpack=True, **kwargs)
    return coords
