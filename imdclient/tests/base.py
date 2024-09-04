from imdclient.IMDClient import IMDClient
import pytest
from pathlib import Path
import os
import signal
import subprocess
import time
import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal
from numpy.testing import (
    assert_array_almost_equal,
    assert_almost_equal,
    assert_allclose,
)

import logging

logger = logging.getLogger("imdclient.IMDClient")


class IMDv3IntegrationTest:

    @pytest.fixture()
    def run_sim_and_wait(self, command, match_string):
        p = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True,
            text=True,
            bufsize=0,
            preexec_fn=os.setsid,
        )
        t = time.time()
        for stdout_line in iter(p.stdout.readline, ""):
            logger.debug(f"stdout: {stdout_line}")
            if match_string in stdout_line:
                break
            if time.time() - t > 10:
                raise TimeoutError("Timeout waiting for match string")

        logger.debug("Match string found")
        yield p
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)

    @pytest.fixture()
    def simulation(self, tmp_path, run_sim_and_wait):
        old_cwd = Path.cwd()
        os.chdir(tmp_path)
        yield
        os.chdir(old_cwd)

    @pytest.fixture()
    def client(self, simulation, universe):
        client = IMDClient("localhost", 8888, universe.trajectory.n_atoms)
        yield client
        client.stop()

    def test_compare_imd_to_true_traj(self, universe, client):
        assert len(universe.trajectory) == 1001
        imdf = client.get_imdframe()
        assert_allclose(imdf.time, universe.trajectory[1].time)
        assert_allclose(imdf.step, universe.trajectory[1].data["step"])
        assert_allclose(imdf.box, universe.trajectory[1].triclinic_dimensions)

        assert_allclose(
            imdf.positions, universe.trajectory[1].positions, rtol=1e-05
        )
        assert_allclose(imdf.forces, universe.trajectory[1].forces, rtol=1e-05)
        assert_allclose(
            imdf.velocities, universe.trajectory[1].velocities, rtol=1
        )
