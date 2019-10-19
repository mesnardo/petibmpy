"""Tests for module `petibmpy.grid`."""

import copy
import numpy
import pathlib
import unittest

import petibmpy


class GridIOTestCase(unittest.TestCase):
    """Tests related to the I/O grid."""

    def setUp(self):
        """Setup."""
        self.x = numpy.sort(numpy.random.rand(5))
        self.y = numpy.sort(numpy.random.rand(10))
        self.z = numpy.sort(numpy.random.rand(20))

    def test_read_write_grid_hdf5(self):
        """Test I/O functions for HDF5 format."""
        filepath = pathlib.Path('grid.h5')
        for dim in [2, 3]:
            coords = [self.x, self.y, self.z][:dim]
            petibmpy.write_grid_hdf5(filepath, 'name', *coords)
            coords2 = petibmpy.read_grid_hdf5(filepath, 'name')
            self.assertEqual(len(coords2), len(coords))
            for i in range(dim):
                self.assertTrue(numpy.allclose(coords2[i], coords[i]))
        filepath.unlink()


class SegmentTestCase(unittest.TestCase):
    """Tests the `Segment` class."""

    def setUp(self):
        """Setup."""
        pass

    def test_create(self):
        """Create a segment given a configuration."""
        # Check uniform segment.
        start, end, width, r, num = 0.0, 1.0, 0.1, 1.0, 11
        config = dict(start=start, end=end, width=width)
        segment = petibmpy.Segment(config=config)
        self.assertEqual(segment.start, start)
        self.assertEqual(segment.end, end)
        self.assertEqual(segment.r, r)
        x = segment.asarray()
        x_true = numpy.linspace(start, end, num=num)
        self.assertTrue(numpy.allclose(x, x_true))
        # Check stretched segment.
        start, end, width, r = 0.0, 1.0, 0.1, 1.01
        config = dict(start=start, end=end, width=width, stretchRatio=r)
        segment = petibmpy.Segment(config=config)
        self.assertEqual(segment.start, start)
        self.assertEqual(segment.end, end)
        self.assertEqual(segment.r, r)
        x = segment.asarray()
        self.assertTrue(abs(x[1] - x[0]) <= width)
        self.assertEqual(abs(x[2] - x[1]) / abs(x[1] - x[0]), segment.r)
        # Check reversed stretched segment.
        start, end, width, r = 0.0, 1.0, 0.1, 1.01
        config = dict(start=start, end=end, width=width, stretchRatio=r,
                      reverse=True)
        segment = petibmpy.Segment(config=config)
        self.assertEqual(segment.start, start)
        self.assertEqual(segment.end, end)
        self.assertAlmostEqual(segment.r, 1 / r, places=12)
        x = segment.asarray()
        self.assertTrue(abs(x[-1] - x[-2]) <= width)
        self.assertAlmostEqual(abs(x[-2] - x[-3]) / abs(x[-1] - x[-2]), r,
                               places=12)


class GridLineTestCase(unittest.TestCase):
    """Tests the `GridLine` class."""

    def setUp(self):
        """Setup."""
        subconfig1 = dict(end=-2.0, width=0.1, stretchRatio=1.01, reverse=True)
        subconfig2 = dict(end=2.0, width=0.1)
        subconfig3 = dict(end=10.0, width=0.1, stretchRatio=1.02)
        subDomains = [subconfig1, subconfig2, subconfig3]
        self.config = dict(direction='x', start=-10.0, subDomains=subDomains)

    def test_create_gridline(self):
        """Create a grid line given a configuration."""
        gridline = petibmpy.GridLine(config=self.config)
        x = gridline.asarray()
        self.assertEqual(x[0], self.config['start'])
        self.assertEqual(x[-1], self.config['subDomains'][-1]['end'])


class CartesianGridTestCase(unittest.TestCase):
    """Tests the `CartesianGrid` class."""

    def setUp(self):
        """Setup."""
        subconfig1 = dict(end=-2.0, width=0.1, stretchRatio=1.01, reverse=True)
        subconfig2 = dict(end=2.0, width=0.1)
        subconfig3 = dict(end=20.0, width=0.1, stretchRatio=1.02)
        config_x = dict(direction='x', start=-10.0,
                        subDomains=[subconfig1, subconfig2, subconfig3])
        config_y = dict(direction='y', start=-15.0,
                        subDomains=[subconfig1, subconfig2, subconfig3])
        config_z = dict(direction='z', start=-2.0, subDomains=[subconfig2])
        self.config = [config_x, config_y, config_z]

    def test_grid_create(self):
        """Create a grid given a configuration."""
        grid = petibmpy.CartesianGrid(config=self.config)
        self.assertEqual(len(grid.gridlines), 3)
        for i, x in enumerate(grid.get_gridlines()):
            config = self.config[i]
            self.assertAlmostEqual(x[0], config['start'], places=12)
            self.assertAlmostEqual(x[-1], config['subDomains'][-1]['end'],
                                   places=12)

    def test_write_yaml(self):
        """Test the method to create the YAML configuration file."""
        grid = petibmpy.CartesianGrid(config=self.config)
        filepath = pathlib.Path('mesh.yaml')
        grid.write_yaml(filepath, ndigits=12)
        with open(filepath, 'r') as infile:
            lines = infile.readlines()
        datadir = pathlib.Path(__file__).absolute().parent / 'data'
        filepath2 = datadir / 'mesh.yaml'
        with open(filepath2, 'r') as infile:
            lines2 = infile.readlines()
        self.assertEqual(len(lines), len(lines2))
        for line, line2 in zip(lines, lines2):
            self.assertEqual(line, line2)
        filepath.unlink()

    def test_write_hdf5(self):
        """Test the method to write the mesh in HDF5 format."""
        grid = petibmpy.CartesianGrid(config=self.config)
        filepath = pathlib.Path('grid.h5')
        grid.write_hdf5(filepath)
        coords = petibmpy.read_grid_hdf5(filepath, 'vertex')
        for x, x2 in zip(grid.get_gridlines(), coords):
            self.assertTrue(numpy.allclose(x, x2))
        filepath.unlink()

    def test_num_cells(self):
        """Test the method to get the number of cells."""
        start, end, width = -10.0, 10.0, 0.1
        num = int(abs(end - start) / width)
        subconfig = dict(start=start, subDomains=[dict(end=end, width=width)])
        config = []
        for dim, direction in zip([1, 2, 3], ['x', 'y', 'z']):
            subconfig['direction'] = direction
            config.append(copy.deepcopy(subconfig))
            grid = petibmpy.CartesianGrid(config=config)
            n_cells = grid.get_number_cells()
            self.assertEqual(n_cells, num**dim)
