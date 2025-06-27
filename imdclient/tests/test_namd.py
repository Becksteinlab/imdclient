import logging
from pathlib import Path
import re

import pytest
from numpy.testing import (
    assert_allclose,
)
import MDAnalysis as mda

from .base import IMDv3IntegrationTest, assert_allclose_with_logging
from .datafiles import (
    NAMD_TOPOL,
    NAMD_CONF_NST_1,
    NAMD_CONF_NST_8,
    NAMD_PARAMS,
    NAMD_PSF,
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

    @pytest.fixture()
    def container_name(self):
        return "ghcr.io/becksteinlab/streaming-namd-docker:main-common-cpu"

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
    def dt(self, inp):
        pattern = re.compile(r"^\s*timestep\s*(\S+)")
        with open(inp, "r") as file:
            for line in file:
                match = pattern.match(line)
                if match:
                    # NAMD timestep is in femtoseconds, convert to picoseconds
                    # as IMDv3 expects dt in ps, 1fs = 0.001ps
                    return float(match.group(1)) * 0.001
        raise ValueError(f"No dt found in {inp}")

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
    def test_compare_imd_to_true_traj(self, imd_u, true_u, first_frame, dt):
        for i in range(first_frame, len(true_u.trajectory)):

            assert_allclose(
                true_u.trajectory[i].time,
                imd_u.trajectory[i - first_frame].time,
                atol=1e-03,
            )

            assert_allclose(
                dt,
                imd_u.trajectory[i - first_frame].dt,
                atol=1e-03,
            )
            # step in DCDReader is frame index, not integration step
            # don't compare step

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

    # Since NAMD does not write velocities, forces to the DCD file, we need to do so seperately by extracting that info from their respective DCD files
    # Compare velocities
    def test_compare_imd_to_true_traj_vel(self, imd_u, true_u_vel, first_frame):
        for i in range(first_frame, len(true_u_vel.trajectory)):

            assert_allclose_with_logging(
                # Unit conversion
                true_u_vel.trajectory[i].positions * 20.45482706,
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
