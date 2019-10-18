"""Tests for module `petibmpy.body`."""

import numpy
import pathlib
import unittest

import petibmpy


class BodyTestCase(unittest.TestCase):
    """Tests for body functions."""

    def setUp(self):
        """Setup."""
        self.datadir = pathlib.Path(__file__).absolute().parent / 'data'
        self.num = 101
        self.x = numpy.random.rand(self.num)
        self.y = numpy.random.rand(self.num)
        self.z = numpy.random.rand(self.num)

    def test_read_body_2d(self):
        """Test function `read_body` for 2D configuration."""
        filepath = self.datadir / 'body2d.txt'
        x, y = petibmpy.read_body(filepath, skiprows=1)
        with open(filepath, 'r') as infile:
            num = int(infile.readline())
        self.assertEqual(x.size, num)
        self.assertEqual(y.size, num)

    def test_read_body_3d(self):
        """Test function `read_body` for 3D configuration."""
        filepath = self.datadir / 'body3d.txt'
        x, y, z = petibmpy.read_body(filepath, skiprows=1)
        with open(filepath, 'r') as infile:
            num = int(infile.readline())
        self.assertEqual(x.size, num)
        self.assertEqual(y.size, num)
        self.assertEqual(z.size, num)

    def test_read_write_body(self):
        """Test functions `read_body` and `write_body`."""
        filepath = self.datadir / 'body.txt'
        coords = (self.x, self.y, self.z)
        for d in [2, 3]:
            petibmpy.write_body(filepath, *coords[:d])
            coords2 = petibmpy.read_body(filepath, skiprows=1)
            self.assertEqual(len(coords2), d)
            for i, coord in enumerate(coords2):
                self.assertTrue(numpy.allclose(coord, coords[i]))
        filepath.unlink()
