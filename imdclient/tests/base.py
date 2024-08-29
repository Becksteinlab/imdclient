import pytest
from pathlib import Path
import os
from .utils import run_sim_and_wait
import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal

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
    def store_imd_traj(self, universe):
        timesteps = []
        for ts in universe.trajectory:
            timesteps.append(ts.copy())
        return timesteps

    @pytest.fixture()
    def true_universe(self, topology, true_trajectory, universe_kwargs):
        return mda.Universe(
            topology,
            true_trajectory,
            **universe_kwargs,
        )

    def test_compare_imd_to_true_traj(self, store_imd_traj, true_universe):
        for i in range(len(true_universe.trajectory)):
            assert_timestep_almost_equal(
                true_universe.trajectory[i], store_imd_traj[i], decimal=3
            )
