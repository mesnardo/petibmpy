"""Module to create/read/write a PetIBM grid."""

import functools
import h5py
import math
from matplotlib import pyplot
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
        self.dim = 0  # number of dimensions
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
            self.dim += 1
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

    def plot_gridlines(self, **kwargs):
        if self.dim == 2:
            return self.plot_gridlines_2d(**kwargs)
        elif self.dim == 3:
            return self.plot_gridlines_3d(**kwargs)
        else:
            raise ValueError(f'dim = {self.dim} not supported!')

    def plot_gridlines_2d(self, figsize=(6.0, 6.0), color='black',
                          xlabel='x', ylabel='y',
                          xrange=(0, None, 1), yrange=(0, None, 1),
                          xlim=(-numpy.infty, numpy.infty),
                          ylim=(-numpy.infty, numpy.infty)):
        """Create a Matplotlib figure with gridlines.

        Parameters
        ----------
        figsize : (float, float), optional
            Width and height of the figure in inches; default is (6, 6).
        color : str, optional
            Color of the gridlines; default is black.
        xlabel : str, optional
            Label along the x axis; default is 'x'.
        ylabel : str, optional
            Label along the y axis; default is 'y'.
        xrange : (int, int, int), optional
            Index range (min, max, stride) to consider for x gridlines;
            default is to consider all stations (0, None, 1).
        yrange : (int, int, int), optional
            Index range (min, max, stride) to consider for y gridlines;
            default is to consider all stations (0, None, 1).
        xlim : (float, float), optional
            Limits of the domain in the x direction to plot;
            default is to plot the entire domain.
        ylim : (float, float), optional
            Limits of the domain in the y direction to plot;
            default is to plot the entire domain.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib Figure.
        matplotlib.axes.Axes
            Matplotlib Axes object.

        """
        assert self.dim == 2
        x, y = self.get_gridlines()
        fig, ax = pyplot.subplots(figsize=figsize)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        self._plot_gridlines_2d(ax, x, y, color=color,
                                xrange=xrange, yrange=yrange,
                                xlim=xlim, ylim=ylim)
        return fig, ax

    def plot_gridlines_3d(self, figsize=(12.0, 6.0), color='black',
                          xlabel='x', ylabel='y', zlabel='z',
                          xrange=(0, None, 1), yrange=(0, None, 1),
                          zrange=(0, None, 1),
                          xlim=(-numpy.infty, numpy.infty),
                          ylim=(-numpy.infty, numpy.infty),
                          zlim=(-numpy.infty, numpy.infty)):
        """Create a Matplotlib figure with gridlines.

        Parameters
        ----------
        figsize : (float, float), optional
            Width and height of the figure in inches; default is (12, 6).
        color : str, optional
            Color of the gridlines; default is black.
        xlabel : str, optional
            Label along the x axis; default is 'x'.
        ylabel : str, optional
            Label along the y axis; default is 'y'.
        zlabel : str, optional
            Label along the z axis; default is 'z'.
        xrange : (int, int, int), optional
            Index range (min, max, stride) to consider for x gridlines;
            default is to consider all stations (0, None, 1).
        yrange : (int, int, int), optional
            Index range (min, max, stride) to consider for y gridlines;
            default is to consider all stations (0, None, 1).
        zrange : (int, int, int), optional
            Index range (min, max, stride) to consider for z gridlines;
            default is to consider all stations (0, None, 1).
        xlim : (float, float), optional
            Limits of the domain in the x direction to plot;
            default is to plot the entire domain.
        ylim : (float, float), optional
            Limits of the domain in the y direction to plot;
            default is to plot the entire domain.
        zlim : (float, float), optional
            Limits of the domain in the z direction to plot;
            default is to plot the entire domain.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib Figure.
        array of matplotlib.axes.Axes
            Array of Matplotlib Axes objects.

        """
        assert self.dim == 3
        x, y, z = self.get_gridlines()
        fig, ax = pyplot.subplots(ncols=3, figsize=figsize)
        ax[0].set_xlabel(xlabel)
        ax[0].set_ylabel(ylabel)
        self._plot_gridlines_2d(ax[0], x, y, color=color,
                                xrange=xrange, yrange=yrange,
                                xlim=xlim, ylim=ylim)
        ax[1].set_xlabel(xlabel)
        ax[1].set_ylabel(zlabel)
        self._plot_gridlines_2d(ax[1], x, z, color=color,
                                xrange=xrange, yrange=zrange,
                                xlim=xlim, ylim=zlim)
        ax[2].set_xlabel(zlabel)
        ax[2].set_ylabel(ylabel)
        self._plot_gridlines_2d(ax[2], z, y, color=color,
                                xrange=zrange, yrange=yrange,
                                xlim=zlim, ylim=ylim)
        return fig, ax

    def _plot_gridlines_2d(self, ax, x, y, color='black',
                           xrange=(0, None, 1), yrange=(0, None, 1),
                           xlim=(-numpy.infty, numpy.infty),
                           ylim=(-numpy.infty, numpy.infty)):
        x = x[xrange[0]:xrange[1]:xrange[2]]
        y = y[yrange[0]:yrange[1]:yrange[2]]
        mask, = numpy.where((x >= xlim[0]) & (x <= xlim[1]))
        x = x[mask]
        mask, = numpy.where((y >= ylim[0]) & (y <= ylim[1]))
        y = y[mask]
        ax.vlines(x, ymin=y[0], ymax=y[-1], color=color)
        ax.hlines(y, xmin=x[0], xmax=x[-1], color=color)
        ax.axis('scaled', adjustable='box')
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(y[0], y[-1])

    def print_info(self):
        """Print some information about the cell widths.

        The method prints the minimum and maximum cell widths
        along each direction, as well as max/min ratio across
        directions.

        """
        gridlines = self.get_gridlines()
        x, y = gridlines[0], gridlines[1]
        dx, dy = x[1:] - x[:-1], y[1:] - y[:-1]
        dx_min, dx_max = numpy.min(dx), numpy.max(dx)
        dy_min, dy_max = numpy.min(dy), numpy.max(dy)
        print(f'dx: min={dx_min}, max={dx_max}')
        print(f'dy: min={dy_min}, max={dy_max}')
        print(f'dx_max / dy_min = {dx_max / dy_min}')
        print(f'dy_max / dx_min = {dy_max / dx_min}')
        if self.dim == 3:
            z = gridlines[2]
            dz = z[1:] - z[:-1]
            dz_min, dz_max = numpy.min(dz), numpy.max(dz)
            print(f'dz: min={dz_min}, max={dz_max}')
            print(f'dx_max / dz_min = {dx_max / dz_min}')
            print(f'dz_max / dx_min = {dz_max / dx_min}')
            print(f'dy_max / dz_min = {dy_max / dz_min}')
            print(f'dz_max / dy_min = {dz_max / dy_min}')


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
            if self._split_needed(node):
                node1, node2 = self._split_uniform_and_stretch(node)
                self.segments.append(Segment(config=node1))
                self.segments.append(Segment(config=node2))
            else:
                self.segments.append(Segment(config=node))
            start = node['end']
        self.size = self.get_size()

    def get_size(self):
        """Return the number of vertices in the gridline."""
        size = self.asarray().size
        return size

    def asarray(self, tol=1e-12):
        """Return the gridline as a 1D NumPy array of floats."""
        vertices = numpy.concatenate(list(s.asarray() for s in self.segments))
        d = numpy.append(True, numpy.diff(vertices))
        return vertices[d > tol]

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

    def _split_needed(self, config):
        """Check if need to split a configuration into uniform and stretched.

        We only to split the configuration is the last width is bigger than
        the target maximum width.

        Parameters
        ----------
        config : dict
            Configuration of the segment to split.

        Returns
        -------
        bool
            True is splitting is needed.

        """
        if 'max_width' not in config.keys():
            return False
        start, end = config['start'], config['end']  # segment limits
        length = abs(end - start)  # segment length
        ratio = config['stretchRatio']  # stretching ratio
        first_width = config['width']  # first width of stretched portion
        max_width = config['max_width']  # largest width of stretched portion
        # Compute number of elements in the geometric series.
        n = round(math.log(1 + length / first_width * (ratio - 1)) /
                  math.log(ratio))
        # Compute the width of the elements.
        widths = numpy.full(n, ratio)
        widths[0] = first_width
        widths = numpy.cumprod(widths)
        # Check if last width is bigger than target maximum width.
        return widths[-1] > max_width

    def _split_uniform_and_stretch(self, config):
        """Split configuration of a stretched portion.

        The configuration is split into a stretch portion and a uniform portion
        with a cell width equal to the maximum cell width provided.

        Parameters
        ----------
        config : dict
            Configuration of the segment to split.

        Returns
        -------
        dict, dict
            Configurations for the stretched and uniform sub-segments.

        """
        start, end = config['start'], config['end']  # segment limits
        length = abs(end - start)  # segment length
        ratio = config['stretchRatio']  # stretching ratio
        first_width = config['width']  # first width of stretched portion
        max_width = config['max_width']  # largest width of stretched portion

        # Find estimation of the position of the junction
        # between the uniform and stretched portion.
        n = 0
        st, width = 0.0, first_width
        while width < max_width:
            width = first_width * ratio**n
            ed = st + width
            st = ed
            n += 1

        # Compute number of cells in uniform portion.
        n = math.ceil(abs(length - ed) / max_width)

        # Create the configuration for the uniform and stretched portions.
        reverse = config.get('reverse', False)
        if reverse:
            junction = start + n * max_width
            config1 = dict(start=start, end=junction, width=max_width)
            config2 = dict(start=junction, end=end, width=first_width,
                           stretchRatio=ratio, reverse=True)
        else:
            junction = end - n * max_width
            config1 = dict(start=start, end=junction, width=first_width,
                           stretchRatio=ratio)
            config2 = dict(start=junction, end=end, width=max_width)

        return config1, config2


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
            assert abs(n - n_float) < 1e-6, "Length should be multiple of width"
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
