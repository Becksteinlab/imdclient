import pytest
from pathlib import Path
import os
from .utils import run_sim_and_wait
import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal
from numpy.testing import assert_array_almost_equal, assert_almost_equal

import logging

logger = logging.getLogger("imdclient.IMDClient")


class IMDv3IntegrationTest:

    @pytest.fixture()
    def simulation(self, tmp_path, command, match_string):
        old_cwd = Path.cwd()
        os.chdir(tmp_path)
        p = run_sim_and_wait(command, match_string)
        yield p
        os.chdir(old_cwd)
        # NOTE: fixme
        p.kill()
        logger.debug("Killed simulation process")
        p.wait()

    @pytest.fixture()
    def universe(self, simulation, universe_kwargs, topology):
        u = mda.Universe(
            topology,
            "localhost:8888",
            **universe_kwargs,
        )
        return u

    @pytest.fixture()
    def true_universe(self, topology, true_trajectory, universe_kwargs):
        return mda.Universe(
            topology,
            true_trajectory,
            **universe_kwargs,
        )

    def test_compare_imd_to_true_traj(self, universe, true_universe):
        assert len(true_universe.trajectory) == 1001
        i = 0
        for i, ts_imd in enumerate(universe.trajectory):
            assert_timestep_almost_equal(ts_imd, true_universe.trajectory[i])

            i += 1
        assert i == 1000
