"""Module for PetIBM probes."""

import collections
import h5py
import numpy
import yaml

from .misc import _setup_yaml, _Sequence


class _ProbeBase(object):
    """Base class for a probe."""

    _type = 'BASE'

    def __init__(self, name, field,
                 viewer='hdf5', path=None,
                 n_monitor=None, n_sum=None, t_start=None, t_end=None):
        """Initialize a base probe.

        Parameters
        ----------
        name : str
            Name of the probe
        field : str
            Name of the field variable to monitor
        viewer : str, optional
            Type of viewer, by default 'hdf5'
        path : pathlib.Path or str, optional
            Path of the output file, by default None
        n_monitor : int, optional
            Monitoring frequency, by default None
        n_sum : int, optional
            Number of time steps to average, by default None
        t_start : float, optional
            Starting time of monitoring, by default None
        t_end : float, optional
            Ending time of monitoring, by default None

        """
        self.name = name
        self.field = field
        self.set_viewer(viewer=viewer, path=path)
        self.n_monitor = n_monitor
        self.n_sum = n_sum
        self.t_start, self.t_end = t_start, t_end

    def __repr__(self):
        """Return the string representation.

        Returns
        -------
        str
            The string representation.

        """
        repr_str = ''
        kwargs = ['n_monitor', 'n_sum', 't_start', 't_end']
        for attr in kwargs:
            if getattr(self, attr) is not None:
                repr_str += f'- {attr}: {getattr(self, attr)}\n'
        return (f'Probe({self.name}):\n' +
                f'- field: {self.field}\n' +
                f'- viewer: {self.viewer}\n' +
                f'- path: {self.path}\n' +
                repr_str)

    def set_viewer(self, viewer='hdf5', path=None):
        """Set the output viewer type and path.

        The path is relative to the PetIBM output directory.

        Parameters
        ----------
        viewer : str, optional
            Type of viewer, choices are 'hdf5' or 'ascii', by default 'hdf5'
        path : pathlib.Path or str, optional
            Path of the output file, by default None

        """
        assert viewer in ['hdf5', 'ascii']
        self.viewer = viewer
        self.path = path

    def _get_yaml_node(self):
        node = collections.OrderedDict({})
        node['name'] = self.name
        node['type'] = self._type
        node['field'] = self.field
        node['viewer'] = self.viewer
        node['path'] = self.path
        for attr in ['n_sum', 'n_monitor', 't_start', 't_end']:
            if getattr(self, attr) is not None:
                node[attr] = getattr(self, attr)
        return node


class ProbeVolume(_ProbeBase):
    """Class for a volume probe (monitoring solution in sub-volume)."""

    _type = 'VOLUME'

    def __init__(self, name, field, box=None,
                 adjust_box=False, grid=None, **kwargs):
        """Initialize a volume probe.

        Parameters
        ----------
        name : str
            Name of the probe
        field : str
            Name of the field variable to monitor
        box : list or numpy.ndarray, optional
            Limits of the box, by default None
        adjust_box : bool, optional
            Adjust the box given a grid, by default False
        grid : list or numpy.ndarray, optional
            The grid of the field, by default None
        **kwargs: dict, optional
            Optional arguments passed to the base class constructor

        """
        super(ProbeVolume, self).__init__(name, field, **kwargs)
        self.box = box
        if box is not None and adjust_box and grid is not None:
            self.adjust_box(grid, box=box)

    def __repr__(self):
        """Return the string representation.

        Returns
        -------
        str
            The string representation

        """
        repr_super = super().__repr__()
        return (repr_super +
                f'- type: {self._type}\n' +
                f'- box: {self.box}\n')

    @classmethod
    def _check_type(cls, ptype):
        """Check if probe type matches class type.

        Parameters
        ----------
        ptype : str
            Type of the probe

        Returns
        -------
        bool
            True if type is 'VOLUME'

        """
        return ptype == ProbeVolume._type

    def _get_yaml_node(self, ndigits=6, **kwargs):
        node = super(ProbeVolume, self)._get_yaml_node(**kwargs)
        assert self.box is not None
        box_node = collections.OrderedDict({})
        dim = len(self.box)
        for i, d in enumerate(['x', 'y', 'z'][:dim]):
            # Round the limits to given number of digits.
            limits = [round(float(l), ndigits=ndigits) for l in self.box[i]]
            box_node[d] = _Sequence(limits)
        node['box'] = box_node
        return node

    def adjust_box(self, grid, box=None):
        """Adjust the box so that limits lie between two grid points.

        Parameters
        ----------
        grid : list or numpy.ndarray
            The grid of the field to minotor
        box : list or numpy.ndarray, optional
            Estimated limits of the box, by default None

        """
        if box is None:
            box = self.box
        self.box = []
        for x, v in zip(grid, box):
            x_l, x_r = v
            idx_l = numpy.where(x < x_l)[0][-1]
            idx_r = numpy.where(x > x_r)[0][0]
            x_l = 0.5 * (x[idx_l] + x[idx_l + 1])
            x_r = 0.5 * (x[idx_r - 1] + x[idx_r])
            self.box.append([x_l, x_r])

    def read_hdf5(self, filepath, time, ndigits=6):
        """Read the probe from a HDF5 file at a given time.

        Parameters
        ----------
        filepath : pathlib.Path or str
            Path of file with the solution of the probe
        time : float
            Time value
        ndigits : int, optional
            Number of digits to round the time value, by default 6

        Returns
        -------
        tuple
            The mesh grid of the probe
        numpy.ndarray
            The probe values

        """
        f = h5py.File(str(filepath), 'r')
        mesh_group = f['mesh']
        dim = len(mesh_group)
        mesh = mesh_group['x'][:], mesh_group['y'][:]
        if dim == 3:
            mesh = mesh + (mesh_group['z'][:],)
        mesh_sizes = [line.size for line in mesh]
        field_group = f[self.field]
        time_str = f'{time:0.{ndigits}f}'
        values = field_group[time_str][:]
        values = values.reshape(mesh_sizes[::-1])
        f.close()
        return mesh, values


class ProbePoint(_ProbeBase):
    """Class to monitor a field at a single point."""

    _type = 'POINT'

    def __init__(self, name, field, loc=None, **kwargs):
        """Initialize a point probe.

        Parameters
        ----------
        name : str
            Name of the probe
        field : str
            Name of the field to monitor
        loc : list or numpy.ndarray, optional
            Coordinates of the point to monitor, by default None
        **kwargs: dict, optional
            Optional arguments passed to the base class constructor

        """
        super(ProbePoint, self).__init__(name, field, **kwargs)
        self.set_loc(loc)

    def __repr__(self):
        """Return the string representation.

        Returns
        -------
        str
            The string representation

        """
        repr_super = super().__repr__()
        return (repr_super +
                f'- type: {self._type}\n' +
                f'- loc: {self.loc}\n')

    @classmethod
    def _check_type(cls, ptype):
        """Check if probe type matches class type.

        Parameters
        ----------
        ptype : str
            Type of the probe

        Returns
        -------
        bool
            True if type is 'POINT'

        """
        return ptype == ProbePoint._type

    def set_loc(self, loc):
        """Set the coordinates of the point to monitor.

        Parameters
        ----------
        loc : list or numpy.ndarray
            Coordinates of the point

        """
        self.loc = loc

    def _get_yaml_node(self, **kwargs):
        node = super(ProbePoint, self)._get_yaml_node(**kwargs)
        if self.loc is not None:
            node['loc'] = _Sequence(self.loc)
        return node


def Probe(ptype, *args, **kwargs):
    """Create a probe.

    Parameters
    ----------
    ptype : str
        Type of the probe, choices are 'VOLUME' or 'POINT'

    Returns
    -------
    ProbeVolume or ProbePoint
        The probe

    Raises
    ------
    ValueError
        Type is neither 'VOLUME' nor 'POINT'

    """
    for subclass in _ProbeBase.__subclasses__():
        if subclass._check_type(ptype):
            return subclass(*args, **kwargs)
    raise ValueError('Parameter `ptype` not in [`VOLUME`, `POINT`]!')


def probes_yaml_dump(probes, filepath, mode='w'):
    """Save the probes configuration in a YAML file.

    Parameters
    ----------
    probes : list
        The list of probes
    filepath : pathlib.Path or str
        Path of the YAML file
    mode : str, optional
        Mode to open file, choices are 'w' or 'a', by default 'w'

    """
    assert mode in ['w', 'a']
    if not hasattr(probes, '__iter__'):
        probes = [probes]
    _setup_yaml()
    probes_yaml = [probe._get_yaml_node() for probe in probes]
    with open(filepath, mode) as outfile:
        yaml.dump({'probes': probes_yaml}, outfile, default_flow_style=False)
