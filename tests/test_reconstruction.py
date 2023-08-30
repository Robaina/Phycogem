#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for module reconstruction.py
"""

import tempfile
import unittest
from pathlib import Path

from phycogem.reconstruction import *

this_file_dir = Path(__file__).parent


class TestGEM(unittest.TestCase):
    def test_remove_shuttle_reactions(self):
        self.assertEqual()


if __name__ == "__main__":
    unittest.main()
