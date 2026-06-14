import logging
from pathlib import Path
import re

import pytest

from .base import IMDv2IntegrationTest, IMDv3IntegrationTest
from .datafiles import (
    GROMACS_GRO,
    GROMACS_TOP,
    GROMACS_MDP_V3_NST_1,
    GROMACS_MDP_V3_NST_8,
    GROMACS_MDP_V2_NST_1,
    GROMACS_MDP_V2_NST_8,
)

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("gromacs_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class IMDGromacsTest:

    @pytest.fixture()
    def setup_command(self, mdp):
        return f"gmx grompp -f {Path(mdp).name} -c {Path(GROMACS_GRO).name} -p {Path(GROMACS_TOP).name} -o topol.tpr"

    @pytest.fixture()
    def simulation_command(self):
        return f"gmx mdrun -s topol.tpr -o ci.trr -imdport 8888 -imdwait"

    @pytest.fixture()
    def input_files(self, mdp):
        return [GROMACS_TOP, mdp, GROMACS_GRO]

    @pytest.fixture()
    def traj(self):
        return "ci.trr"

    @pytest.fixture()
    def topol(self):
        return Path(GROMACS_GRO).name

    @pytest.fixture()
    def dt(self, mdp):
        pattern = re.compile(r"^\s*dt\s*=\s*(\S+)")
        with open(mdp, "r") as file:
            for line in file:
                match = pattern.match(line)
                if match:
                    return float(match.group(1))
        raise ValueError(f"No dt found in {mdp}")

    # @pytest.fixture()
    # def match_string(self):
    #     return "IMD: Will wait until I have a connection and IMD_GO orders."

    @pytest.fixture()
    def first_frame(self):
        return 0


class TestIMDv3Gromacs(IMDGromacsTest, IMDv3IntegrationTest):

    @pytest.fixture(params=[GROMACS_MDP_V3_NST_1, GROMACS_MDP_V3_NST_8])
    def mdp(self, request):
        return request.param


class TestIMDv2Gromacs(IMDGromacsTest, IMDv2IntegrationTest):

    @pytest.fixture(params=[GROMACS_MDP_V2_NST_1, GROMACS_MDP_V2_NST_8])
    def mdp(self, request):
        return request.param

    @pytest.fixture()
    def post_simulation_command(self):
        # Match GROMACS IMD v2 mol-shifted whole coords using pbc whole.
        return (
            "echo '0' | gmx trjconv -f ci.trr -s topol.tpr -o ci_whole.trr -pbc whole && "
            "mv ci_whole.trr ci.trr"
        )
