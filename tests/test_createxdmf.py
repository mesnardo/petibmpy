"""Tests for module `petibmpy.body`."""

import filecmp
import numpy
import pathlib
import unittest

import petibmpy


class CreateXDMFTestCase(unittest.TestCase):
    """Test functions in the module `createxdmf`."""

    def setUp(self):
        """Setup."""
        datadir = pathlib.Path(__file__).absolute().parent / 'data'
        self.gridpath_2d = datadir / 'grid2d.h5'
        self.gridpath_3d = datadir / 'grid3d.h5'

    def test_write_xdmf_2d(self):
        """Test the function `write_xdmf` for a 2D configuration."""
        names = ['p', 'u', 'v']
        for name in names:
            datadir = pathlib.Path('none')
            nstart, nt, nsave = 0, 1000, 100
            filepath1 = pathlib.Path(name + '1.xmf')
            petibmpy.write_xdmf(filepath1, datadir, self.gridpath_2d, name,
                                nstart=nstart, nt=nt, nsave=nsave)
            filepath2 = pathlib.Path(name + '2.xmf')
            states = list(range(nstart, nt + 1, nsave))
            petibmpy.write_xdmf(filepath2, datadir, self.gridpath_2d, name,
                                states=states)
            self.assertTrue(filecmp.cmp(filepath1, filepath2))
            filepath1.unlink()
            filepath2.unlink()

    def test_write_xdmf_3d(self):
        """Test the function `write_xdmf` for a 3D configuration."""
        names = ['p', 'u', 'v', 'w']
        for name in names:
            datadir = pathlib.Path('none')
            nstart, nt, nsave = 0, 1000, 100
            filepath1 = pathlib.Path(name + '1.xmf')
            petibmpy.write_xdmf(filepath1, datadir, self.gridpath_3d, name,
                                nstart=nstart, nt=nt, nsave=nsave)
            filepath2 = pathlib.Path(name + '2.xmf')
            states = list(range(nstart, nt + 1, nsave))
            petibmpy.write_xdmf(filepath2, datadir, self.gridpath_3d, name,
                                states=states)
            self.assertTrue(filecmp.cmp(filepath1, filepath2))
            filepath1.unlink()
            filepath2.unlink()

    def test_write_xdmf_multi_2d(self):
        """Test the function `write_xdmf_multi` for a 2D configuration."""
        datadir = pathlib.Path('none')
        config = dict(grid=self.gridpath_2d,
                      data={'p': datadir, 'u': datadir, 'v': datadir})
        nstart, nt, nsave = 0, 1000, 100
        filepath1 = pathlib.Path('1.xmf')
        petibmpy.write_xdmf_multi(filepath1, config,
                                  nstart=nstart, nt=nt, nsave=nsave)
        filepath2 = pathlib.Path('2.xmf')
        states = list(range(nstart, nt + 1, nsave))
        petibmpy.write_xdmf_multi(filepath2, config, states=states)
        self.assertTrue(filecmp.cmp(filepath1, filepath2))
        filepath1.unlink()
        filepath2.unlink()

    def test_write_xdmf_multi_3d(self):
        """Test the function `write_xdmf_multi` for a 3D configuration."""
        datadir = pathlib.Path('none')
        config = dict(grid=self.gridpath_3d,
                      data={'p': datadir, 'u': datadir,
                            'v': datadir, 'w': datadir})
        nstart, nt, nsave = 0, 1000, 100
        filepath1 = pathlib.Path('1.xmf')
        petibmpy.write_xdmf_multi(filepath1, config,
                                  nstart=nstart, nt=nt, nsave=nsave)
        filepath2 = pathlib.Path('2.xmf')
        states = list(range(nstart, nt + 1, nsave))
        petibmpy.write_xdmf_multi(filepath2, config, states=states)
        self.assertTrue(filecmp.cmp(filepath1, filepath2))
        filepath1.unlink()
        filepath2.unlink()
