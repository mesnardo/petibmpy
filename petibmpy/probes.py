"""Module for the PetIBM probes."""

import collections
import h5py
import numpy
import yaml

from .misc import setup_yaml, Sequence


class ProbeBase(object):

    def __init__(self, name, field, viewer='hdf5', path=None,
                 n_monitor=None, n_sum=None, t_start=None, t_end=None):
        self.name = name
        self.field = field
        assert viewer in ['hdf5', 'ascii']
        self.viewer = viewer
        self.path = path
        self.n_monitor = n_monitor
        self.n_sum = n_sum
        self.t_start, self.t_end = t_start, t_end

    def _get_yaml_node(self):
        node = collections.OrderedDict({})
        node['name'] = self.name
        node['type'] = self._type
        node['field'] = self.field
        node['viewer'] = self.viewer
        node['relative_path'] = self.path
        for attr in ['n_sum', 'n_monitor', 't_start', 't_end']:
            if getattr(self, attr) is not None:
                node[attr] = getattr(self, attr)
        return node


class ProbeVolume(ProbeBase):
    _type = 'VOLUME'

    def __init__(self, *args, **kwargs):
        super(ProbeVolume, self).__init__(*args, **kwargs)
        self.box = collections.OrderedDict({'x': None, 'y': None, 'z': None})

    @classmethod
    def check_type(cls, ptype):
        return ptype == ProbeVolume._type

    def _get_yaml_node(self):
        node = super(ProbeVolume, self)._get_yaml_node()
        box_node = collections.OrderedDict({})
        dim = len(self.box)
        for d in ['x', 'y', 'z'][:dim]:
            box_node[d] = Sequence(self.box[d])
        node['box'] = box_node
        return node

    def read_hdf5(self, filepath, name, time):
        f = h5py.File(str(filepath), 'r')
        mesh_group = f['mesh']
        dim = len(mesh_group)
        mesh = mesh_group['x'][:], mesh_group['y'][:]
        if dim == 3:
            mesh = mesh + (mesh_group['z'][:],)
        mesh_sizes = [line.size for line in mesh]
        field_group = f[name]
        time_str = '{:0.6f}'.format(round(time, ndigits=6))
        values = field_group[time_str][:]
        values = values.reshape(mesh_sizes[::-1])
        f.close()
        return mesh, values

    def set_box(self, grid, volume, ndigits=6):
        self.dim = len(grid)
        dirs = ['x', 'y', 'z'][:self.dim]
        self.box = collections.OrderedDict({})
        for d, x, v in zip(dirs, grid, volume):
            x_l, x_r = v
            idx_l = numpy.where(x < x_l)[0][-1]
            idx_r = numpy.where(x > x_r)[0][0]
            x_l = 0.5 * (x[idx_l] + x[idx_l + 1])
            x_r = 0.5 * (x[idx_r - 1] + x[idx_r])
            x_l = round(float(x_l), ndigits=ndigits)
            x_r = round(float(x_r), ndigits=ndigits)
            self.box[d] = [x_l, x_r]


class ProbePoint(ProbeBase):
    _type = 'POINT'

    def __init__(self, *args, **kwargs):
        super(ProbePoint, self).__init__(*args, **kwargs)
        self.loc = [None, None, None]

    @classmethod
    def check_type(cls, ptype):
        return ptype == ProbePoint._type

    def _get_yaml_node(self):
        node = super(ProbePoint, self)._get_yaml_node()
        node['loc'] = Sequence(self.loc)
        return node


def Probe(ptype, *args, **kwargs):
    for subclass in ProbeBase.__subclasses__():
        if subclass.check_type(ptype):
            return subclass(*args, **kwargs)
    raise ValueError('Parameter `ptype` not in [`VOLUME`, `POINT`]!')


def probes_yaml_dump(probes, filepath, mode='w'):
    assert mode in ['w', 'a']
    if not hasattr(probes, '__iter__'):
        probes = [probes]
    setup_yaml()
    probes_yaml = [probe._get_yaml_node() for probe in probes]
    with open(filepath, mode) as outfile:
        yaml.dump({'probes': probes_yaml}, outfile, default_flow_style=False)
