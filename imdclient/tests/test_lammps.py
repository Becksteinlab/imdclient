import MDAnalysis as mda
import pytest
import subprocess
from pathlib import Path
import os
import logging
from .utils import run_sim_and_wait, get_free_port
from .base import IMDv3IntegrationTest
from .datafiles import LAMMPS_IN, LAMMPS_TOPOL

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
    def topology(self):
        return LAMMPS_TOPOL

    @pytest.fixture()
    def traj_path(self):
        """Relative path of the trajectory file within `tmp_dir`.
        Allows comparison of large trajectory files without storing them in the repo.
        """
        return Path("trj.dump")

    @pytest.fixture()
    def universe_kwargs(self):
        """Simulation-engine specific keyword arguments for the Universe."""
        return {
            "atom_style": "id type x y z",
            "timeout": 10,
        }
