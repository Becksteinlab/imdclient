"""
IMDClient
"""

# Don't import IMDReader here, eventually it may be moved to a separate package
from .IMDClient import IMDClient
from importlib.metadata import version

from .streamanalysis import AnalysisBase, StackableAnalysis
from MDAnalysis.analysis import base

base.AnalysisBase = AnalysisBase

__version__ = version("imdclient")
