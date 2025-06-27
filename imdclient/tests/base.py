import time
import logging
import shutil

import pytest
import numpy as np
from numpy.testing import (
    assert_allclose,
)
import docker
import MDAnalysis as mda

from .utils import get_free_port
from .minimalreader import MinimalReader

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
    def container_name(self):
        return "ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu"

    @pytest.fixture()
    def setup_command(self):
        return None

    @pytest.fixture()
    def port(self):
        yield get_free_port()

    @pytest.fixture()
    def docker_client(
        self,
        tmp_path,
        input_files,
        setup_command,
        simulation_command,
        port,
        container_name,
    ):
        # In CI, container process needs access to tmp_path
        tmp_path.chmod(0o777)
        docker_client = docker.from_env()
        img = docker_client.images.pull(
            container_name,
        )
        # Copy input files into tmp_path
        for inp in input_files:
            shutil.copy(inp, tmp_path)

        cmdstring = "cd '/tmp'"
        # Run the setup command, if any
        if setup_command is not None:
            # This should be blocking
            cmdstring += " && " + setup_command

        cmdstring += " && " + simulation_command

        # Start the container, mount tmp_path, run simulation
        container = docker_client.containers.run(
            img,
            f"/bin/sh -c '{cmdstring}'",
            detach=True,
            volumes={tmp_path.as_posix(): {"bind": "/tmp", "mode": "rw"}},
            ports={"8888/tcp": port},
            name="sim",
            remove=True,
        )

        # For now, just wait 30 seconds
        # life is too short to figure out how to redirect all stdout from inside
        # a container
        time.sleep(30)

        yield
        try:
            container.stop()
        except docker.errors.NotFound:
            pass

    @pytest.fixture()
    def imd_u(self, docker_client, topol, tmp_path, port):
        n_atoms = mda.Universe(tmp_path / topol).atoms.n_atoms
        u = MinimalReader(
            f"imd://localhost:{port}", n_atoms=n_atoms, process_stream=True
        )
        yield u

    @pytest.fixture()
    def true_u(self, topol, traj, imd_u, tmp_path):
        u = mda.Universe(
            (tmp_path / topol),
            (tmp_path / traj),
        )
        yield u

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

            assert_allclose(
                true_u.trajectory[i].data["step"],
                imd_u.trajectory[i - first_frame].step,
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

            assert_allclose_with_logging(
                true_u.trajectory[i].velocities,
                imd_u.trajectory[i - first_frame].velocities,
                atol=1e-03,
            )

            assert_allclose_with_logging(
                true_u.trajectory[i].forces,
                imd_u.trajectory[i - first_frame].forces,
                atol=1e-03,
            )

    def test_continue_after_disconnect(
        self, docker_client, topol, tmp_path, port
    ):
        n_atoms = mda.Universe(
            tmp_path / topol,
            # Make sure LAMMPS topol can be read
            # Does nothing if not LAMMPS
            atom_style="id type x y z",
        ).atoms.n_atoms
        u = MinimalReader(
            f"imd://localhost:{port}",
            n_atoms=n_atoms,
            continue_after_disconnect=True,
        )
        # Though we disconnect here, the simulation should continue
        u.close()
        # Wait for the simulation to finish running
        time.sleep(45)

        # Now, attempt to reconnect- should fail,
        # since the simulation should have continued
        with pytest.raises(IOError):
            n_atoms = mda.Universe(
                tmp_path / topol,
                atom_style="id type x y z",
            ).atoms.n_atoms
            u = MinimalReader(f"imd://localhost:{port}", n_atoms=n_atoms)

    def test_wait_after_disconnect(self, docker_client, topol, tmp_path, port):
        n_atoms = mda.Universe(
            tmp_path / topol,
            # Make sure LAMMPS topol can be read
            # Does nothing if not LAMMPS
            atom_style="id type x y z",
        ).atoms.n_atoms
        u = MinimalReader(
            f"imd://localhost:{port}",
            n_atoms=n_atoms,
            continue_after_disconnect=False,
        )
        u.close()
        # Give the simulation engine
        # enough time to finish running (though it shouldn't)
        time.sleep(45)

        n_atoms = mda.Universe(
            tmp_path / topol,
            atom_style="id type x y z",
        ).atoms.n_atoms
        u = MinimalReader(f"imd://localhost:{port}", n_atoms=n_atoms)
