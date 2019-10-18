"""Run the test suite."""

import sys
import unittest


tests = ['test_rotate', 'test_forces', 'test_body', 'test_extrude']

suite = unittest.TestSuite()

for test in tests:
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(test))

unittest.TextTestRunner().run(suite).wasSuccessful()
