"""
IMDClient
"""

# Don't import IMDReader here, eventually it may be moved to a separate package
from .IMDClient import IMDClient
from importlib.metadata import version

__version__ = version("imdclient")
