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
