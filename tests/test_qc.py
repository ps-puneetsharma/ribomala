#!/usr/bin/env python3

"""
Tests for the ribosome profiling quality control functionality.
"""

import shutil
from pathlib import Path

import pysam
import polars as pl
import pytest

# Import the module to test
from ribomala import qc

