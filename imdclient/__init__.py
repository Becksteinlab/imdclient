"""
IMDReader
"""

from .IMDClient import IMDClient
from .IMDREADER import IMDReader
from importlib.metadata import version


__version__ = version("imdclient")
