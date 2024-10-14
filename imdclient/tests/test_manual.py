from imdclient.IMDClient import IMDClient
import pytest
import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal
from numpy.testing import (
    assert_array_almost_equal,
    assert_almost_equal,
    assert_allclose,
)
import numpy as np
from .base import assert_allclose_with_logging

import logging

logger = logging.getLogger("imdclient.IMDClient")


class TestIMDv3Manual:
    """
    Tool for running IMDv3 integration tests via the command line.

    To use, start the simulation, wait for it to be ready for an IMD connection,
    and then run this command relative to the root of the cloned respository:

    pytest -s imdclient/tests/test_manual.py \
        --topol_arg <path/to/topology> \
        --traj_arg <path/to/trajectory> \
        --first_frame_arg <first traj frame to compare to IMD>

    Where the topology is the same topology as the IMD system, the trajectory is the trajectory
    to compare to IMD data read from the socket, and the first frame is the first frame of the 
    trajectory which should be compared to IMD data read from the socket (0 for GROMACS and NAMD, 1 for LAMMPS)
    """

    @pytest.fixture()
    def universe(self, topol_arg, traj_arg):
        return mda.Universe(topol_arg, traj_arg)

    @pytest.fixture()
    def client(self, universe):
        client = IMDClient("localhost", 8888, universe.trajectory.n_atoms)
        yield client
        client.stop()

    def test_compare_imd_to_true_traj(self, universe, client, first_frame_arg):
        imdsinfo = client.get_imdsessioninfo()

        for ts in universe.trajectory[first_frame_arg:]:
            imdf = client.get_imdframe()
            if imdsinfo.time:
                assert_allclose(imdf.time, ts.time, atol=1e-03)
                assert_allclose(imdf.step, ts.data["step"])
            if imdsinfo.box:
                assert_allclose_with_logging(
                    imdf.box,
                    ts.triclinic_dimensions,
                    atol=1e-03,
                )
            if imdsinfo.positions:
                assert_allclose_with_logging(
                    imdf.positions, ts.positions, atol=1e-03
                )
            if imdsinfo.velocities:
                assert_allclose_with_logging(
                    imdf.velocities, ts.velocities, atol=1e-03
                )
            if imdsinfo.forces:
                assert_allclose_with_logging(
                    imdf.forces, ts.forces, atol=1e-03
                )
