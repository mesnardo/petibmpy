"""Fake tests."""

import unittest

import petibmpy


class FakeTest(unittest.TestCase):
    """Fake tests."""

    def setUP(self):
        """Setup."""
        pass

    def test_print_version(self):
        """Print the package version."""
        print(petibmpy.__version__)


if __name__ == '__main__':
    unittest.main()
