from imdclient.IMDClient import IMDClient
from imdclient.IMDREADER import IMDReader
import pytest
from pathlib import Path
import os
import signal
import subprocess
import time
from numpy.testing import (
    assert_allclose,
)
import numpy as np
import docker
import logging
import shutil
import MDAnalysis as mda

logger = logging.getLogger("imdclient.IMDClient")


def assert_allclose_with_logging(a, b, rtol=1e-07, atol=0, equal_nan=False):
    """
    Custom function to compare two arrays element-wise, similar to np.testing.assert_allclose,
    but logs all non-matching values.

    Parameters:
    a, b : array_like
        Input arrays to compare.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.
    equal_nan : bool
        Whether to compare NaNs as equal.
    """
    # Convert inputs to numpy arrays
    a = np.asarray(a)
    b = np.asarray(b)

    # Compute the absolute difference
    diff = np.abs(a - b)

    # Check if values are within tolerance
    not_close = diff > (atol + rtol * np.abs(b))

    # Check if there are any NaNs and handle them if necessary
    if equal_nan:
        nan_mask = np.isnan(a) & np.isnan(b)
        not_close &= ~nan_mask

    # Log all the values that are not close
    if np.any(not_close):
        print("The following values do not match within tolerance:")
        for idx in np.argwhere(not_close):
            logger.debug(
                f"a[{tuple(idx)}]: {a[tuple(idx)]}, b[{tuple(idx)}]: {b[tuple(idx)]}, diff: {diff[tuple(idx)]}"
            )
        # Optionally raise an error after logging if you want it to behave like assert
        raise AssertionError("Arrays are not almost equal.")
    else:
        print("All values are within tolerance.")


class IMDv3IntegrationTest:

    @pytest.fixture()
    def setup_command(self):
        return None

    @pytest.fixture()
    def docker_client(
        self,
        tmp_path,
        input_files,
        setup_command,
        simulation_command,
        match_string,
    ):
        docker_client = docker.from_env()
        docker_client.images.pull(
            "ghcr.io/becksteinlab/streaming-md-docker:main"
        )
        # Copy input files into tmp_path
        for inp in input_files:
            shutil.copy(inp, tmp_path)
        # Start the container, mount tmp_path
        container = docker_client.containers.run(
            "ghcr.io/becksteinlab/streaming-md-docker:main",
            volumes={tmp_path: {"bind": "/tmp", "mode": "rw"}},
            detach=True,
        )
        # Run the setup command, if any
        if setup_command is not None:
            # This should be blocking
            container.exec_run(setup_command, workdir="/tmp")

        # Start the simulation
        cmd = docker_client.api.exec_create(
            container.id, simulation_command, workdir="/tmp"
        )
        exec_output = docker_client.api.exec_start(exec_id=cmd["Id"])
        # Wait for the match string
        wait_time = 0
        t = time.time()
        matched = False
        while wait_time < 60:
            try:
                # Open and read the log file directly from the host-side mounted path
                with open(tmp_path / "output.log", "r") as log_file:
                    log_content = log_file.read()

                logger.debug(
                    f"{tmp_path / "output.log"} content: {log_content.strip()}"
                )
                # Check if the match_string exists in the log content
                if match_string in log_content:
                    matched = True
                    break

            except FileNotFoundError:
                wait_time = time.time() - t
                time.sleep(1)

            time.sleep(1)

        if not matched:
            raise RuntimeError(
                f"Simulation failed to reach IMD Session readiness"
            )

        yield

        container.stop()

    @pytest.fixture()
    def imd_u(self, docker_client, topol, tmp_path):
        u = mda.Universe((tmp_path / topol), "imd://localhost:8888")
        with mda.Writer((tmp_path / "imd.trr"), u.trajectory.n_atoms) as w:
            for ts in u.trajectory:
                w.write(u.atoms)
        yield mda.Universe((tmp_path / topol), (tmp_path / "imd.trr"))

    def test_compare_imd_to_true_traj(self, imd_u, true_u, first_frame):
        for i in range(first_frame, len(true_u.trajectory)):
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
