"""
IMDReader
"""

from .IMDClient import IMDClient
from .IMDREADER import IMDReader
from importlib.metadata import version

from .streamanalysis import AnalysisBase, StackableAnalysis
from MDAnalysis.analysis import base

base.AnalysisBase = AnalysisBase

__version__ = version("imdclient")
