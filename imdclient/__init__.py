"""
IMDReader
"""

from .IMDClient import IMDClient
from importlib.metadata import version


__version__ = version("imdclient")
