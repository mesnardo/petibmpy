"""Module to create/read/write a PetIBM grid."""

import functools
import h5py
import math
import numpy
import operator
import yaml


class CartesianGrid():
    """Contain information about a structured Cartesian grid."""

    def __init__(self, config=None):
        """Initialize the grid.

        Parameters
        ----------
        config : dictionary (optional)
            Configuration of the grid to create; default: None.

        """
        self.gridlines = {'x': GridLine(), 'y': GridLine(), 'z': GridLine()}
        self.n = 0  # number of cells
        if config is not None:
            self.create(config)

    def __repr__(self, ndigits=6):
        """Representation of the grid.

        Parameters
        ----------
        ndigits : integer (optional)
            Number of digits to represent floats; default: 6.

        """
        sub_repr = ',\n'.join((line.__repr__(ndigits=ndigits)
                               for line in self.gridlines.values()))
        return ('Grid(cells={}, gridlines=[\n{}])'.format(self.n, sub_repr))

    def create(self, config):
        """Create the grid.

        Parameters
        ----------
        config : dictionary
            Configuration of the grid.

        """
        gridlines = {}
        for node in config:
            direction = node['direction']
            assert direction in self.gridlines.keys()
            gridlines[direction] = GridLine(config=node)
        self.gridlines = gridlines
        self.n = self.get_number_cells()

    def get_number_cells(self):
        """Return the number of cells in the grid."""
        n_gridlines = [line.size - 1 for line in self.gridlines.values()]
        n = functools.reduce(operator.mul, n_gridlines, 1)
        return n

    def get_gridlines(self):
        """Return the gridlines as a list of 1D NumPy arrays of floats."""
        gridlines = []
        for direction in ['x', 'y', 'z']:
            if direction in self.gridlines.keys():
                gridlines.append(self.gridlines[direction].asarray())
        return gridlines

    def write_hdf5(self, filepath):
        """Save the grid into HDF5 file.

        Parameters
        ----------
        filepath : pathlib.Path or string
            Path of the HDF5 file to write into.

        """
        gridlines = self.get_gridlines()
        write_grid_hdf5(filepath, 'vertex', *gridlines)

    def write_yaml(self, filepath, ndigits=6):
        """Write the YAML configuration node for PetIBM.

        Parameters
        ----------
        filepath : pathlib.Path or string
            Path of the YAML file to write into.
        ndigits : integer (optional)
            Number of digits to represent floats; default: 6.

        """
        n_gridlines = [line.size - 1 for line in self.gridlines.values()]
        nodes = [gridline.yaml_node(ndigits=ndigits)
                 for gridline in self.gridlines.values()]
        with open(filepath, 'w') as outfile:
            outfile.write('# ' +
                          'x'.join(map(str, n_gridlines)) +
                          ' ({})\n'.format(self.n))
            outfile.write(yaml.dump({'mesh': nodes}, default_flow_style=False))


class GridLine():
    """Contain information about a gridline of a structured Cartesian grid."""

    def __init__(self, config=None):
        """Initialize the gridline.

        Parameters
        ----------
        config : dictionary (optional)
            Configuration of the gridline to create; default: None.

        """
        self.direction = None
        self.segments = []
        self.size = 0  # number of vertices
        if config is not None:
            self.create(config)

    def __repr__(self, ndigits=6):
        """Representation of the gridline.

        Parameters
        ----------
        ndigits : integer (optional)
            Number of digits to represent floats; default: 6.

        """
        sub_repr = ',\n'.join((s.__repr__(ndigits=ndigits)
                               for s in self.segments))
        return ('Gridline(direction={}, size={}, segments=[\n{}])'
                .format(self.direction, self.size, sub_repr))

    def create(self, config):
        """Create the gridline.

        Parameters
        ----------
        config : dictionary
            Configuration of the gridline.

        """
        self.direction = config['direction']
        start = config['start']
        for node in config['subDomains']:
            node['start'] = start
            self.segments.append(Segment(config=node))
            start = node['end']
        self.size = self.get_size()

    def get_size(self):
        """Return the number of vertices in the gridline."""
        size = self.asarray().size
        return size

    def asarray(self):
        """Return the gridline as a 1D NumPy array of floats."""
        vertices = numpy.concatenate(list(s.asarray() for s in self.segments))
        return numpy.unique(vertices)

    def yaml_node(self, ndigits=6):
        """Return the YAML configuration node for PetIBM.

        Parameters
        ----------
        ndigits : integer (optional)
            Number of digits to represent floats; default: 6.

        Returns
        -------
        node : dictionary
            Configuration node for the gridline.

        """
        start = self.segments[0].asarray()[0]
        node = {'direction': self.direction,
                'start': float(round(start, ndigits=ndigits)),
                'subDomains': [s.yaml_node(ndigits=ndigits)
                               for s in self.segments]}
        return node


class Segment():
    """Contain information about a segment of a gridline."""

    def __init__(self, config=None):
        """Initialize the segment.

        Parameters
        ----------
        config : dictionary (optional)
            Configuration of the segment to create; default: None.

        """
        self.vertices = []
        self.start, self.end = None, None
        self.n = 0  # number of vertices
        self.r = 1.0  # stretching ratio
        if config is not None:
            self.create(config)

    def __repr__(self, ndigits=6):
        """Representation of the segment.

        Parameters
        ----------
        ndigits : integer (optional)
            Number of digits to represent floats; default: 6.

        """
        return ('Segment(start={}, end={}, n={}, r={})'
                .format(round(self.start, ndigits=ndigits),
                        round(self.end, ndigits=ndigits),
                        self.n,
                        round(self.r, ndigits=ndigits)))

    def create(self, config):
        """Create the segment.

        Parameters
        ----------
        config : dictionary
            Configuration of the segment.

        """
        start, end = config['start'], config['end']
        length = abs(end - start)
        width = config['width']
        r = config.get('stretchRatio', 1.0)
        reverse = config.get('reverse', False)

        if abs(r - 1) < 1e-6:
            n_float = length / width
            n = int(round(n_float))
            assert abs(n - n_float) < 1e-6, "Length should be multple of width"
            self.vertices = numpy.linspace(start, end, num=n + 1)
            self.r = 1.0
        else:
            n_float = math.log(1 + length / width * (r - 1)) / math.log(r)
            n = int(round(n_float))
            width = length * (r - 1) / (r**n - 1)  # re-compute first width
            widths = [width * r**k for k in range(n)]
            cumsum = numpy.cumsum(widths)
            if reverse:
                self.vertices = numpy.concatenate(([end], end - cumsum))[::-1]
                self.r = abs((self.vertices[-1] - self.vertices[-2]) /
                             (self.vertices[-2] - self.vertices[-3]))
            else:
                self.vertices = numpy.concatenate(([start], start + cumsum))
                self.r = abs((self.vertices[2] - self.vertices[1]) /
                             (self.vertices[1] - self.vertices[0]))
        self.start, self.end = start, end
        self.n = self.vertices.size

    def asarray(self):
        """Return the segment as a 1D NumPy array of floats."""
        return numpy.array(self.vertices)

    def yaml_node(self, ndigits=6):
        """Return the YAML configuration node for PetIBM.

        Parameters
        ----------
        ndigits : integer (optional)
            Number of digits to represent floats; default: 6.

        Returns
        -------
        node : dictionary
            Configuration node for the segment.

        """
        node = {'end': float(round(self.end, ndigits=ndigits)),
                'cells': len(self.vertices) - 1,
                'stretchRatio': float(round(self.r, ndigits=ndigits))}
        return node


def read_grid_hdf5(filepath, name):
    """Read a grid from HDF5 file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the HDF5 file.
    name : string
        Name of the grid.

    Returns
    -------
    x : numpy.ndarray
        The x-coordinates along a gridline in the x-direction.
    y : numpy.ndarray
        The y-coordinates along a gridline in the y-direction.
    z : numpy.ndarray
        The z-coordinates along a gridline in the z-direction.

    """
    f = h5py.File(str(filepath), 'r')
    dim = len(f[name])
    x, y, z = f[name]['x'][:], f[name]['y'][:], None
    if dim == 3:
        z = f[name]['z'][:]
    f.close()
    if z is None or len(z) == 1:
        return x, y
    return x, y, z


def write_grid_hdf5(filepath, name, *grid):
    """Write a grid to a HDF5 file.

    Parameters
    ----------
    filepath : string or pathlib.Path object
        Path of the HDF5 file.
    name : string
        Name of the grid.
    grid : tuple of numpy.ndarray objects
        The gridline coordinates as 1D arrays of floats.

    """
    labels = ('x', 'y', 'z')
    f = h5py.File(str(filepath), 'w')
    group = f.create_group(name)
    for i, gridline in enumerate(grid):
        group.create_dataset(labels[i], data=gridline)
    f.close()
    return
