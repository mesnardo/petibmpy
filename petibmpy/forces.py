"""Module with functions to process forces.."""

import numpy


def read_forces(filepath):
    """Read PetIBM forces from given file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the file to read.

    Returns
    -------
    data : tuple of 4 numpy.ndarray objects
        Time followed by the forces in the x, y, and z directions.

    """
    with open(filepath, 'r') as infile:
        data = numpy.loadtxt(infile, unpack=True)
    return data


def get_force_coefficients(*forces, coeff=1.0):
    """Convert forces to force coefficients.

    Parameters
    ----------
    forces : tuple of numpy.ndarray objects
        The forces.
    coeff : float (optional)
        The scaling coefficient; default: 1.0.

    Returns
    -------
    force_coeffs : tuple of numpy.ndarray objects
        The force coefficients.

    """
    force_coeffs = (coeff * f for f in forces)
    return force_coeffs


def get_time_averaged_values(t, *forces, limits=(-numpy.infty, numpy.infty)):
    """Compute the time-averaged values.

    Parameters
    ----------
    t : numpy.ndarray object
        The time values.
    forces : tuple of numpy.ndarray objects
        The forces (or force coefficients).
    limits : tuple of 2 floats (optional)
        Time limits used to compute the mean; default: (-inf, +inf).

    Returns
    -------
    means : tuple of floats
        The time-averaged values.

    """
    mask = (t >= limits[0]) & (t <= limits[1])
    means = (numpy.mean(f[mask]) for f in forces)
    return means
