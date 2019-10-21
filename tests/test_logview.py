"""Tests for module `petibmpy.logview`."""

from matplotlib import pyplot
import pathlib
import unittest

import petibmpy


class PETScLogViewTestCase(unittest.TestCase):
    """Tests for functions about PETSc log view."""

    def setUp(self):
        """Setup."""
        datadir = pathlib.Path(__file__).absolute().parent / 'data'
        self.filepath = datadir / 'view.log'

    def test_parse_log_view(self):
        """Test function `parse_log_view`."""
        log_view = petibmpy.PETScLogView(self.filepath)
        # Check wall-time clock and resident set size.
        self.assertEqual(log_view.walltime, 1425.0)
        self.assertIsNone(log_view.res)
        # Check name of events.
        event_names = ['Main Stage', 'initialize', 'rhsVelocity',
                       'solveVelocity', 'rhsPoisson', 'solvePoisson',
                       'update', 'write', 'monitor']
        for name in log_view.events.keys():
            self.assertIn(name, event_names)
        # Check event about solving the Poisson system.
        event = log_view.events['solvePoisson']
        self.assertEqual(event['index'], 5)
        self.assertEqual(event['wall-time (s)'], 1066.3)
        self.assertEqual(event['wall-time (%)'], 74.8)

    def test_plot_events_breakdown(self):
        """Test function `plot_events_breakdown`."""
        log_view = petibmpy.PETScLogView(filepath=self.filepath)
        events = {'Run 1': log_view.events, 'Run 2': log_view.events}
        fig, ax = pyplot.subplots()
        petibmpy.plot_events_breakdown(ax, events)
