from imdclient.IMDREADER import IMDReader
import pytest
import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal
from numpy.testing import (
    assert_allclose,
)
import numpy as np
from .base import assert_allclose_with_logging
from pathlib import Path

import logging

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("manual_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3Manual:
    """
    Tool for running IMDv3 integration tests via the command line.

    To use, start the simulation, wait for it to be ready for an IMD connection,
    and then run this command relative to the root of the cloned respository:

    pytest -s imdclient/tests/test_manual.py \
        --topol_path_arg <path/to/topology> \
        --traj_path_arg <path/to/trajectory> \
        --first_frame_arg <first traj frame to compare to IMD>

    Where the topology is the same topology as the IMD system, the trajectory is the path where
    the trajectory of the running simulation is being written, and the first frame is the first frame of the
    trajectory which should be compared to IMD data read from the socket (0 for GROMACS and NAMD, 1 for LAMMPS)
    """

    @pytest.fixture()
    def true_u(self, imd_u, topol_path_arg, traj_path_arg):
        return mda.Universe(topol_path_arg, traj_path_arg)

    @pytest.fixture()
    def imd_u(self, topol_path_arg, tmp_path):
        tmp_u = mda.Universe(topol_path_arg, "imd://localhost:8888")
        with mda.Writer(
            f"{tmp_path.as_posix()}/imd_test_traj.trr", tmp_u.atoms.n_atoms
        ) as w:
            for ts in tmp_u.trajectory:
                w.write(tmp_u.atoms)
        imd_u = mda.Universe(
            topol_path_arg, f"{tmp_path.as_posix()}/imd_test_traj.trr"
        )
        yield imd_u

    def test_compare_imd_to_true_traj(self, true_u, imd_u, first_frame_arg):

        for i in range(first_frame_arg, len(true_u.trajectory)):
            assert_allclose(
                true_u.trajectory[i].time, imd_u.trajectory[i].time, atol=1e-03
            )
            assert_allclose(
                true_u.trajectory[i].data["step"],
                imd_u.trajectory[i].data["step"],
            )
            if true_u.trajectory[i].dimensions is not None:
                assert_allclose_with_logging(
                    true_u.trajectory[i].dimensions,
                    imd_u.trajectory[i].dimensions,
                    atol=1e-03,
                )
            if true_u.trajectory[i].has_positions:
                assert_allclose_with_logging(
                    true_u.trajectory[i].positions,
                    imd_u.trajectory[i].positions,
                    atol=1e-03,
                )
            if true_u.trajectory[i].has_velocities:
                assert_allclose_with_logging(
                    true_u.trajectory[i].velocities,
                    imd_u.trajectory[i].velocities,
                    atol=1e-03,
                )
            if true_u.trajectory[i].has_forces:
                assert_allclose_with_logging(
                    true_u.trajectory[i].forces,
                    imd_u.trajectory[i].forces,
                    atol=1e-03,
                )
