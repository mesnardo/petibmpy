"""Tests for module `petibmpy.regularize`."""

import numpy
import unittest

import petibmpy


class RegularizeTestCase(unittest.TestCase):
    """Tests for regularization function."""

    def setUp(self):
        """Setup."""
        self.num = 20
        theta = numpy.linspace(0.0, 2 * numpy.pi, num=self.num + 1)[:-1]
        self.x, self.y = numpy.cos(theta), numpy.sin(theta)
        self.ds = numpy.sqrt((self.x[1] - self.x[0])**2 +
                             (self.y[1] - self.y[0])**2)

    def test_regularize2d(self):
        """Test function to regularize a 2D curve."""
        # If no optional arguments are provided, curve should be unchanged.
        x, y = petibmpy.regularize2d(self.x, self.y)
        self.assertTrue(numpy.allclose(x, self.x))
        self.assertTrue(numpy.allclose(y, self.y))
        # If curve is already regularized and same number of points is given,
        # curve should be unchanged.
        x, y = petibmpy.regularize2d(self.x, self.y, N=self.num)
        self.assertTrue(numpy.allclose(x, self.x))
        self.assertTrue(numpy.allclose(y, self.y))
        # If curve is already regularized and same spacing is given,
        # curve should be unchanged.
        x, y = petibmpy.regularize2d(self.x, self.y, ds=self.ds)
        self.assertTrue(numpy.allclose(x, self.x))
        self.assertTrue(numpy.allclose(y, self.y))
        # Check size of regularized curve given number of divisions.
        for N in [10, 30]:
            x, y = petibmpy.regularize2d(self.x, self.y, N=N)
            self.assertEqual(x.size, N)
            self.assertEqual(y.size, N)
        # Check spacing of regularized curve given spacing.
        for ds in [self.ds / 2, 2 * self.ds]:
            x, y = petibmpy.regularize2d(self.x, self.y, ds=ds)
            ds2 = numpy.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
            self.assertAlmostEqual(ds2, ds, places=6)
