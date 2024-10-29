import MDAnalysis as mda
import pytest
import logging
from .base import IMDv3IntegrationTest
from .datafiles import NAMD_TOPOL, NAMD_CONF, NAMD_PARAMS, NAMD_PSF
from pathlib import Path

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("namd_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3NAMD(IMDv3IntegrationTest):

    @pytest.fixture()
    def simulation_command(self):
        return f"namd3 {Path(NAMD_CONF).name}"

    @pytest.fixture()
    def input_files(self):
        return [NAMD_TOPOL, NAMD_CONF, NAMD_PARAMS, NAMD_PSF]

    @pytest.fixture()
    def topol(self):
        return Path(NAMD_TOPOL).name

    @pytest.fixture()
    def traj(self):
        return "alanin.dcd"

    @pytest.fixture()
    def comp_dt(self):
        return False

    @pytest.fixture()
    def comp_step(self):
        return False

    @pytest.fixture()
    def comp_time(self):
        return False

    # @pytest.fixture()
    # def match_string(self):
    #     return "INTERACTIVE MD AWAITING CONNECTION"

    @pytest.fixture()
    def first_frame(self):
        return 0


class TestIMDv3NAMDVelocities(IMDv3IntegrationTest):

    @pytest.fixture()
    def simulation_command(self):
        return f"namd3 {Path(NAMD_CONF).name}"

    @pytest.fixture()
    def input_files(self):
        return [NAMD_TOPOL, NAMD_CONF, NAMD_PARAMS, NAMD_PSF]

    @pytest.fixture()
    def topol(self):
        return Path(NAMD_TOPOL).name

    @pytest.fixture()
    def traj(self):
        return "alanin.velocity.dcd"

    @pytest.fixture()
    def comp_dt(self):
        return False

    @pytest.fixture()
    def comp_step(self):
        return False

    @pytest.fixture()
    def comp_time(self):
        return False

    # @pytest.fixture()
    # def match_string(self):
    #     return "INTERACTIVE MD AWAITING CONNECTION"

    @pytest.fixture()
    def first_frame(self):
        return 0


class TestIMDv3NAMDForces(IMDv3IntegrationTest):

    @pytest.fixture()
    def simulation_command(self):
        return f"namd3 {Path(NAMD_CONF).name}"

    @pytest.fixture()
    def input_files(self):
        return [NAMD_TOPOL, NAMD_CONF, NAMD_PARAMS, NAMD_PSF]

    @pytest.fixture()
    def topol(self):
        return Path(NAMD_TOPOL).name

    @pytest.fixture()
    def traj(self):
        return "alanin.force.dcd"

    @pytest.fixture()
    def comp_dt(self):
        return False

    @pytest.fixture()
    def comp_step(self):
        return False

    @pytest.fixture()
    def comp_time(self):
        return False

    # @pytest.fixture()
    # def match_string(self):
    #     return "INTERACTIVE MD AWAITING CONNECTION"

    @pytest.fixture()
    def first_frame(self):
        return 0
