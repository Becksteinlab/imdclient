import MDAnalysis as mda
import pytest
import logging
from .base import IMDv3IntegrationTest, assert_allclose_with_logging
from .datafiles import (
    NAMD_TOPOL,
    NAMD_CONF_NST_1,
    NAMD_CONF_NST_8,
    NAMD_PARAMS,
    NAMD_PSF,
)
from pathlib import Path
from numpy.testing import (
    assert_allclose,
)

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("namd_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3NAMD(IMDv3IntegrationTest):

    @pytest.fixture(params=[NAMD_CONF_NST_1, NAMD_CONF_NST_8])
    def inp(self, request):
        return request.param

    @pytest.fixture()
    def simulation_command(self, inp):
        return f"namd3 {Path(inp).name}"

    @pytest.fixture()
    def input_files(self, inp):
        return [NAMD_TOPOL, inp, NAMD_PARAMS, NAMD_PSF]

    @pytest.fixture()
    def topol(self):
        return Path(NAMD_TOPOL).name

    @pytest.fixture()
    def true_u(self, topol, imd_u, tmp_path):
        u = mda.Universe(
            (tmp_path / topol),
            (tmp_path / "alanin.dcd"),
        )
        yield u

    @pytest.fixture()
    def true_u_vel(self, topol, imd_u, tmp_path):
        u = mda.Universe(
            (tmp_path / topol),
            (tmp_path / "alanin.vel.dcd"),
        )
        yield u

    @pytest.fixture()
    def true_u_force(self, topol, imd_u, tmp_path):
        u = mda.Universe(
            (tmp_path / topol),
            (tmp_path / "alanin.force.dcd"),
        )
        yield u

    # @pytest.fixture()
    # def match_string(self):
    #     return "INTERACTIVE MD AWAITING CONNECTION"

    @pytest.fixture()
    def first_frame(self):
        return 0

    # Compare coords, box, time, dt, step
    def test_compare_imd_to_true_traj(self, imd_u, true_u, first_frame):
        for i in range(first_frame, len(true_u.trajectory)):
            assert_allclose(
                true_u.trajectory[i].time,
                imd_u.trajectory[i - first_frame].time,
                atol=1e-03,
            )
            assert_allclose(
                true_u.trajectory[i].dt,
                imd_u.trajectory[i - first_frame].dt,
                atol=1e-03,
            )
            assert_allclose(
                true_u.trajectory[i].data["step"],
                imd_u.trajectory[i - first_frame].data["step"],
            )
            assert_allclose_with_logging(
                true_u.trajectory[i].dimensions,
                imd_u.trajectory[i - first_frame].dimensions,
                atol=1e-03,
            )
            assert_allclose_with_logging(
                true_u.trajectory[i].positions,
                imd_u.trajectory[i - first_frame].positions,
                atol=1e-03,
            )

    # Compare velocities
    def test_compare_imd_to_true_traj_vel(
        self, imd_u, true_u_vel, first_frame
    ):
        for i in range(first_frame, len(true_u_vel.trajectory)):
            assert_allclose_with_logging(
                true_u_vel.trajectory[i].positions,
                imd_u.trajectory[i - first_frame].velocities,
                atol=1e-03,
            )

    # Compare forces
    def test_compare_imd_to_true_traj_forces(
        self, imd_u, true_u_force, first_frame
    ):
        for i in range(first_frame, len(true_u_force.trajectory)):
            assert_allclose_with_logging(
                true_u_force.trajectory[i].positions,
                imd_u.trajectory[i - first_frame].forces,
                atol=1e-03,
            )
