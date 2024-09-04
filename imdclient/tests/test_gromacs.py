import MDAnalysis as mda
import pytest
from pathlib import Path
import logging

from .utils import get_free_port
from .base import IMDv3IntegrationTest
from .datafiles import GROMACS_TOPOL, GROMACS_TRAJ

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3Gromacs(IMDv3IntegrationTest):

    @pytest.fixture()
    def command(self):
        return f"gmx mdrun -s {GROMACS_TOPOL} -nt 1 -deffnm run -imdport 8888 -imdwait"

    @pytest.fixture()
    def match_string(self):
        return "IMD: Will wait until I have a connection and IMD_GO orders."

    @pytest.fixture()
    def universe(self):
        return mda.Universe(
            GROMACS_TOPOL,
            GROMACS_TRAJ,
        )
