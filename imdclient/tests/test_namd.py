import MDAnalysis as mda
import pytest
import logging
from .base import IMDv3IntegrationTest
from .datafiles import NAMD_TOPOL, NAMD_CONF, NAMD_TRAJ, NAMD_PARAMS, NAMD_PSF

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3NAMD(IMDv3IntegrationTest):

    @pytest.fixture()
    def command(self):
        return (
            f"cp {NAMD_PARAMS} {NAMD_PSF} {NAMD_TOPOL} . && namd3 {NAMD_CONF}"
        )

    @pytest.fixture()
    def match_string(self):
        return "INTERACTIVE MD AWAITING CONNECTION"

    @pytest.fixture()
    def first_frame(self):
        return 0

    @pytest.fixture()
    def universe(self):
        return mda.Universe(
            NAMD_TOPOL,
            NAMD_TRAJ,
        )
