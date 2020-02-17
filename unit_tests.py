# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:42:31 2020
This meant to be used for the tests of units & modules.
@author: vladg
"""

from main_sim import calcCircum
import unittest

class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_calcCircum(self):
        assert calcCircum(0,0) == 0


if __name__ == '__main__':
    unittest.main()