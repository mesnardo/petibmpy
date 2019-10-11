"""Tests for module `petibmpy.forces`."""

import pathlib
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


if __name__ == '__main__':
    unittest.main()
