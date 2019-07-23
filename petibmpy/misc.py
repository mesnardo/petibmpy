"""Collection of miscellaneous functions and classes."""

import collections
import yaml


class Sequence(list):
    """Dummy class to store list/tuple in YAML file in pretty format."""

    pass


def _represent_dictionary_order(self, dict_data):
    """Pretty output of dictionary to YAML file."""
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())


def _represent_limits(self, data):
    """Pretty output of list/tuple to YAML file."""
    return self.represent_sequence('tag:yaml.org,2002:seq', data,
                                   flow_style=True)


def setup_yaml():
    """Configure output format to YAML file."""
    yaml.add_representer(collections.OrderedDict, _represent_dictionary_order)
    yaml.add_representer(Sequence, _represent_limits)
