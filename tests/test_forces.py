"""Tests for module `petibmpy.forces`."""

import numpy
import pathlib
import random
import unittest

import petibmpy


class ForcesTestCase(unittest.TestCase):
    """Tests for forces functions."""

    def setUp(self):
        """Setup."""
        self.datadir = pathlib.Path(__file__).absolute().parent / 'data'

    def test_read_forces(self):
        """Test function `read_forces`."""
        # Tests for files containing forces for single and multiple bodies.
        for n, name in enumerate(['single', 'multiple']):
            dim = [2, 3]  # test 2D and 3D
            for d in dim:
                filepaths = []
                for i in range(2):  # test read from single and multiple files
                    filepath = self.datadir / f'forces{d}d-{name}-{i * 5}.txt'
                    filepaths.append(filepath)
                    res = petibmpy.read_forces(*filepaths)
                    self.assertTrue(len(res) == d * (n + 1) + 1)
                    num = res[0].size
                    self.assertTrue(all(r.size == num for r in res))

    def test_get_force_coefficients(self):
        """Test function `get_force_coefficients`."""
        dim = 3
        # Load forces from file.
        filepath = self.datadir / f'forces{dim}d-single-0.txt'
        _, fx, fy, fz = petibmpy.read_forces(filepath)
        # Convert one direction with unit coefficient.
        cd, = petibmpy.get_force_coefficients(fx)
        self.assertTrue(numpy.allclose(cd, fx))
        # Convert two directions with unit coefficient.
        cd, cl = petibmpy.get_force_coefficients(fx, fy)
        self.assertTrue(numpy.allclose(cd, fx))
        self.assertTrue(numpy.allclose(cl, fy))
        # Convert three directions with unit coefficient.
        cd, cl, cz = petibmpy.get_force_coefficients(fx, fy, fz)
        self.assertTrue(numpy.allclose(cd, fx))
        self.assertTrue(numpy.allclose(cl, fy))
        self.assertTrue(numpy.allclose(cz, fz))
        # Convert three directions with random coefficient.
        coeff = random.uniform(0.1, 5.0)
        cd, cl, cz = petibmpy.get_force_coefficients(fx, fy, fz, coeff=coeff)
        self.assertTrue(numpy.allclose(cd / coeff, fx))
        self.assertTrue(numpy.allclose(cl / coeff, fy))
        self.assertTrue(numpy.allclose(cz / coeff, fz))

    def test_get_time_averaged_values(self):
        """Test function `get_time_averaged_values`."""
        dim = 3
        # Load forces from file.
        filepath1 = self.datadir / f'forces{dim}d-single-0.txt'
        filepath2 = self.datadir / f'forces{dim}d-single-5.txt'
        t, fx, fy, fz = petibmpy.read_forces(filepath1, filepath2)
        # Compute time-averaged value in one direction.
        fx1, = petibmpy.get_time_averaged_values(t, fx)
        self.assertAlmostEqual(fx1, numpy.mean(fx), places=7)
        # Compute time-averaged value for each direction.
        fx1, fy1, fz1 = petibmpy.get_time_averaged_values(t, fx, fy, fz)
        self.assertAlmostEqual(fx1, numpy.mean(fx), places=7)
        self.assertAlmostEqual(fy1, numpy.mean(fy), places=7)
        self.assertAlmostEqual(fz1, numpy.mean(fz), places=7)
        # Compute averaged values providing larger time interval.
        time_limits = (-1.0, 100.0)
        fx1, fy1, fz1 = petibmpy.get_time_averaged_values(t, fx, fy, fz,
                                                          limits=time_limits)
        self.assertAlmostEqual(fx1, numpy.mean(fx), places=7)
        self.assertAlmostEqual(fy1, numpy.mean(fy), places=7)
        self.assertAlmostEqual(fz1, numpy.mean(fz), places=7)
        # Compute averaged values providing restricted time interval.
        time_limits = (5.0, 6.0)
        fx1, fy1, fz1 = petibmpy.get_time_averaged_values(t, fx, fy, fz,
                                                          limits=time_limits)
        self.assertAlmostEqual(fx1, 0.55)
        self.assertAlmostEqual(fy1, 1.55)
        self.assertAlmostEqual(fz1, 2.55)
        # Check if runtime error is raised when wrong time limits are set.
        time_limits = (5.5, 5.5)  # 5.5 is not a saved time value
        with self.assertRaises(RuntimeError):
            fx1, = petibmpy.get_time_averaged_values(t, fx, limits=time_limits)

    def test_get_rms_values(self):
        """Test function `get_rms_values`."""
        a = 2.1
        t = numpy.linspace(0.0, 2 * numpy.pi, num=101)[:-1]
        y1 = numpy.full_like(t, a)
        y2 = a * numpy.sin(t)
        rms1, = petibmpy.get_rms_values(t, y1)
        self.assertAlmostEqual(rms1, a)
        rms1, rms2 = petibmpy.get_rms_values(t, y1, y2)
        self.assertAlmostEqual(rms1, a)
        self.assertAlmostEqual(rms2, a / numpy.sqrt(2))
        rms1, = petibmpy.get_rms_values(t, y1, limits=(0.0, t[t.size // 2]))
        self.assertAlmostEqual(rms1, a)
        limits = (10.0, 10.0)  # 5.5 is not a saved time value
        with self.assertRaises(RuntimeError):
            rms1, = petibmpy.get_rms_values(t, y1, limits=limits)
