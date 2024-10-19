import MDAnalysis as mda
import pytest
import logging
from .base import IMDv3IntegrationTest
from .datafiles import GROMACS_GRO, GROMACS_MDP, GROMACS_TOP
from pathlib import Path

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("gromacs_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3Gromacs(IMDv3IntegrationTest):

    @pytest.fixture()
    def setup_command(self):
        return f"gmx grompp -f {Path(GROMACS_MDP).name} -c {Path(GROMACS_GRO).name} -p {Path(GROMACS_TOP).name} -o topol.tpr"

    @pytest.fixture()
    def simulation_command(self):
        return f"gmx mdrun -s topol.tpr -o ci.trr -imdport 8888 -imdwait"

    @pytest.fixture()
    def input_files(self):
        return [GROMACS_TOP, GROMACS_MDP, GROMACS_GRO]

    @pytest.fixture()
    def traj(self):
        return "ci.trr"

    @pytest.fixture()
    def topol(self):
        return Path(GROMACS_GRO).name

    # @pytest.fixture()
    # def match_string(self):
    #     return "IMD: Will wait until I have a connection and IMD_GO orders."

    @pytest.fixture()
    def first_frame(self):
        return 0
