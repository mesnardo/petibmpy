"""Tests for module `petibmpy.probes`."""


import numpy
import pathlib
import unittest

import petibmpy


class ProbesTestCase(unittest.TestCase):
    """Test the `probes` module."""

    def setUp(self):
        """Setup."""
        loc = (1.1, -2.2, 3.3)
        self.point = petibmpy.Probe('POINT', 'point', 'p', loc=loc)
        box = ((-1.0, 1.1), (-2.0, 2.1), (-3.0, 3.1))
        self.volume = petibmpy.Probe('VOLUME', 'volume', 'p', box=box)

    def test_probes_yaml_dump(self):
        """Test the function `probes_yaml_dump`."""
        filepath = pathlib.Path('probes.yaml')
        petibmpy.probes_yaml_dump([self.point, self.volume], filepath)
        filepath.unlink()
