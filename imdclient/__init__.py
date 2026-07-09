"""
IMDClient
"""

from .IMDClient import IMDClient
from .IMDServer import IMDServer
from importlib.metadata import version

__version__ = version("imdclient")
