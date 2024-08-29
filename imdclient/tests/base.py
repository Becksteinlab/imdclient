import pytest
from pathlib import Path
import os
from .utils import run_sim_and_wait
import MDAnalysis as mda


class IMDv3IntegrationTest:

    @pytest.fixture()
    def simulation(self, tmp_path, command, match_string):
        old_cwd = Path.cwd()
        os.chdir(tmp_path)
        p = run_sim_and_wait(command, match_string)
        yield p
        os.chdir(old_cwd)
        p.kill()
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
    def store_imd_traj(self, universe, topology, traj_path, universe_kwargs):
        timesteps = []
        for ts in universe.trajectory:
            timesteps.append(ts.copy())
        return timesteps

    @pytest.fixture()
    def true_universe(
        self, topology, traj_path, universe_kwargs, store_imd_traj
    ):
        return mda.Universe(
            topology,
            traj_path,
            **universe_kwargs,
        )

    def test_compare_imd_to_true_traj(self, store_imd_traj, true_universe):
        for i, (imd_ts, true_ts) in enumerate(
            zip(store_imd_traj, true_universe.trajectory)
        ):
            assert imd_ts == true_ts
