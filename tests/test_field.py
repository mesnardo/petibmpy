"""Tests for module `petibmpy.field`."""

import numpy
import pathlib
import random
from scipy import interpolate
import unittest

import petibmpy


class FieldTestCase(unittest.TestCase):
    """Test the module `field`."""

    def setUp(self):
        """Setup."""
        self.name = 'p'
        nx, ny, nz = 30, 20, 10
        self.values = numpy.random.rand(nz, ny, nx)
        self.filepath = pathlib.Path(self.name + '.h5')

    def test_read_write_field_hdf5(self):
        """Test the I/O functions."""
        petibmpy.write_field_hdf5(self.filepath, self.name, self.values)
        values = petibmpy.read_field_hdf5(self.filepath, self.name)
        self.assertTrue(numpy.allclose(values, self.values))
        self.filepath.unlink()

    def test_linear_interpolation(self):
        """Test the function to do linear interpolation along first axis."""
        nz = numpy.shape(self.values)[0]
        z = numpy.linspace(0.0, 1.0, num=nz)
        i = -2
        # Exact location.
        values = petibmpy.linear_interpolation(self.values, z, z[i])
        self.assertTrue(numpy.allclose(values, self.values[i]))
        # Mid location.
        zi = 0.5 * (z[i - 1] + z[i])
        values = petibmpy.linear_interpolation(self.values, z, zi)
        self.assertTrue(numpy.allclose(values,
                                       0.5 * (self.values[i - 1] +
                                              self.values[i])))
        # Compare results to scipy.interpolate.interp1d.
        nx = 100
        x = numpy.linspace(0.0, 1.0, num=nx)
        values = numpy.random.rand(nx)
        xi = random.random()
        v1 = petibmpy.linear_interpolation(values, x, xi)
        f = interpolate.interp1d(x, values)
        self.assertAlmostEqual(v1, f([xi])[0], places=12)

    def test_interpolate2d(self):
        """Test the function to interpolate a 2D field onto a grid."""
        _, ny, nx = numpy.shape(self.values)
        x = numpy.linspace(-1.0, 1.0, num=nx)
        y = numpy.linspace(-5.0, 5.0, num=ny)
        values = self.values[0]
        # Interpolate on the same grid.
        values2 = petibmpy.interpolate2d(values, (x, y), (x, y))
        self.assertTrue(numpy.allclose(values2, values))

    def test_interpolate3d(self):
        """Test the function to interpolate a 3D field onto a grid."""
        nz, ny, nx = numpy.shape(self.values)
        x = numpy.linspace(-1.0, 1.0, num=nx)
        y = numpy.linspace(-5.0, 5.0, num=ny)
        z = numpy.linspace(-10.0, 10.0, num=nz)
        # Interpolate on the same grid.
        values2 = petibmpy.interpolate3d(self.values, (x, y, z), (x, y, z))
        self.assertTrue(numpy.allclose(values2, self.values))
