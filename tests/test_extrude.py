"""Tests for module `petibmpy.extrude`."""

import numpy
import random
import unittest

import petibmpy


class ExtrudeTestCase(unittest.TestCase):
    """Tests for the function to extrude 2D bodies."""

    def setUp(self):
        """Setup."""
        self.num = 10
        self.x = numpy.random.rand(self.num)
        self.y = numpy.random.rand(self.num)

    def test_extrude2d(self):
        """Test function `extrude2d`."""
        limits = [-0.5, 0.5]
        n = 5
        # Check fail if neither `n` nor `ds` are provided.
        with self.assertRaises(ValueError):
            petibmpy.extrude2d(self.x, self.y, limits)
        # Check fail if limits provided at too close.
        left_limit = 1.0
        right_limit = left_limit + random.random() * 1e-6
        with self.assertRaises(ValueError):
            petibmpy.extrude2d(self.x, self.y, [left_limit, right_limit], n=n)
        # Check size of 3D body
        x, y, z = petibmpy.extrude2d(self.x, self.y, limits, n=n)
        self.assertEqual(x.size, n * self.x.size)
        self.assertEqual(y.size, n * self.y.size)
        self.assertEqual(z.size, x.size)
        self.assertTrue(numpy.allclose(x[:self.x.size], self.x))
        # Check size of 3D body when given spacing ds.
        ds = abs(limits[1] - limits[0]) / n
        x, y, z = petibmpy.extrude2d(self.x, self.y, limits, ds=ds)
        self.assertEqual(x.size, n * self.x.size)
        self.assertEqual(y.size, n * self.y.size)
        self.assertEqual(z.size, x.size)
        self.assertTrue(numpy.allclose(x[:self.x.size], self.x))
        # Check `force` optional argument.
        x, y, z = petibmpy.extrude2d(self.x, self.y, limits, n=n, force=True)
        self.assertEqual(z[0], limits[0])
        self.assertEqual(z[-1], limits[-1])
