import MDAnalysis as mda
import pytest
import subprocess
from pathlib import Path
import os
import logging

from .utils import get_free_port
from .base import IMDv3IntegrationTest
from .datafiles import LAMMPS_IN, LAMMPS_TOPOL, LAMMPS_TRAJ

# NOTE: removeme
from imdclient.IMDREADER import IMDReader


logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3Lammps(IMDv3IntegrationTest):

    @pytest.fixture()
    def command(self):
        return f"lmp < {LAMMPS_IN}"

    @pytest.fixture()
    def match_string(self):
        return "Waiting for IMD connection on port 8888"

    @pytest.fixture()
    def universe(self):
        return mda.Universe(
            LAMMPS_TOPOL,
            LAMMPS_TRAJ,
            atom_style="id type x y z",
            convert_units=False,
        )
