"""Collection of miscellaneous functions and classes."""

import collections
import functools
import h5py
import inspect
import yaml


def check_not_primary_variables(f):
    """Check if variable names are not primary variables."""
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        forbidden_names = ('vertex', 'u', 'v', 'w', 'p')
        func_args = inspect.getcallargs(f, *args, **kwargs)
        names = func_args.get('names', [])
        common_names = list(set(names) & set(forbidden_names))
        if len(common_names) > 0:
            raise ValueError('Attempting to delete primary variables '
                             f'{common_names}')
        return f(*args, **kwargs)
    return wrapper


@check_not_primary_variables
def delete_datasets_hdf5(filepath, names):
    """Delete datasets from HDF5 file.

    If a name if not a dataset, the function moves to the next name.

    Parameters
    ----------
    filepath : pathlib.Path or str
        Path of the HDF5 file.
    names : list or tuple
        Names of the datasets to delete.

    """
    with h5py.File(filepath, 'a') as f:
        for name in names:
            if name in f.keys():
                del f[name]


class _Sequence(list):
    """Dummy class to store list/tuple in YAML file in pretty format."""

    pass


def _represent_dictionary_order(self, dict_data):
    """Pretty output of dictionary to YAML file."""
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())


def _represent_limits(self, data):
    """Pretty output of list/tuple to YAML file."""
    return self.represent_sequence('tag:yaml.org,2002:seq', data,
                                   flow_style=True)


def _setup_yaml():
    """Configure output format to YAML file."""
    yaml.add_representer(collections.OrderedDict, _represent_dictionary_order)
    yaml.add_representer(_Sequence, _represent_limits)
