"""Tests for module `petibmpy.rotate`."""

import numpy
import unittest

import petibmpy


class RotateTestCase(unittest.TestCase):
    """Tests for rotation functions."""

    def setUp(self):
        """Setup."""
        R = 0.5
        xc, yc = 1.0, 2.0
        theta = numpy.linspace(0.0, 2 * numpy.pi, num=101)[:-1]
        self.x = xc + R * numpy.cos(theta)
        self.y = yc + R * numpy.sin(theta)

    def test_rotate2d(self):
        """Test function `rotate2d`."""
        x, y = petibmpy.rotate2d(self.x, self.y)
        self.assertTrue(numpy.allclose(x, self.x))
        self.assertTrue(numpy.allclose(y, self.y))
        x, y = petibmpy.rotate2d(self.x, self.y, angle=360.0)
        self.assertTrue(numpy.allclose(x, self.x))
        self.assertTrue(numpy.allclose(y, self.y))
        x, y = petibmpy.rotate2d(self.x, self.y,
                                 angle=2 * numpy.pi, mode='rad')
        self.assertTrue(numpy.allclose(x, self.x))
        self.assertTrue(numpy.allclose(y, self.y))
        with self.assertRaises(ValueError):
            petibmpy.rotate2d(x, y, mode='typo')


if __name__ == '__main__':
    unittest.main()
