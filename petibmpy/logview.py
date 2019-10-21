"""Module to parse a PETSc log view file."""

import pathlib
import re


class PETScLogView(object):
    """Parse a PETSc log view file."""

    def __init__(self, filepath=None):
        """Initialize the parser."""
        self.filepath = filepath  # path of the log file
        self.walltime = None  # wall-time clock
        self.res = None  # resident set size
        self.events = {}  # information about logged events
        if self.filepath is not None:
            self.parse_log_view(filepath)

    def parse_log_view(self, filepath):
        """Parse a PETSc log view file."""
        _ = self._read_walltime()
        _ = self._read_resident_set_size()
        _ = self._read_events()

    def _read_walltime(self):
        """Parse and return the wall-time clock in seconds."""
        search = 'Time (sec):'
        with open(self.filepath, 'r') as infile:
            for line in infile.readlines():
                if search in line:
                    self.walltime = float(line.split()[2])
                    break
        if self.walltime is None:
            print('WARNING: Wall-time clock not found.')
        return self.walltime

    def _read_resident_set_size(self, unit='GB'):
        """Parse and return the resident set size."""
        search = 'process memory'
        units = {'B': 0, 'KB': 1, 'MB': 2, 'GB': 3}
        with open(self.filepath, 'r') as infile:
            for line in infile.readlines():
                if search in line:
                    self.res = float(line.split()[7]) / 1024**units[unit]
                    break
        if self.res is None:
            print('WARNING: Resident set size not found.')
        return self.res

    def _read_events(self):
        """Parse information about PETSc events."""
        search = 'Summary of Stages'
        with open(self.filepath, 'r') as infile:
            lines = infile.readlines()
            for i, line in enumerate(lines):
                if search in line:
                    i += 2
                    while not re.match(r'^\s*$', lines[i]):
                        name, data = self._parse_event(lines[i])
                        self.events[name] = data
                        i += 1
                    break
        return self.events

    def _parse_event(self, line):
        """Parse information about an event."""
        info = re.split(':', re.sub(' +', ' ', line))
        name = info[1].strip()
        data = {}
        data['index'] = int(info[0])
        data['wall-time (s)'] = float(info[2].split()[0])
        data['wall-time (%)'] = float(info[2].split()[1][:-1])
        data['FLOPS'] = float(info[2].split()[2])
        return name, data


def plot_events_breakdown(ax, runs,
                          ylabel='wall-time (s)', event_names=None,
                          bar_width=0.5):
    """Add a bar chart of the breakdown of events to an axis."""
    ax.set_ylabel(ylabel)
    index = 1
    for name, events in runs.items():
        if event_names is not None:
            events = dict((name, events[name]) for name in event_names)
        i, offset = 0, 0.0
        for label, event in events.items():
            value = event[ylabel]
            ax.bar(index, value, bottom=offset, label=label,
                   width=bar_width, color=f'C{i}')
            offset += value
            i += 1
        if index == 1:
            ax.legend(loc='upper right', frameon=False)
        index += 1
    ax.set_xticks(list(range(1, index)))
    ax.set_xticklabels(runs.keys())
