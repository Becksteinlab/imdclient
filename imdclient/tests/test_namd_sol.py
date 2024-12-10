from imdclient.IMDREADER import IMDReader
import pytest
import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal
from numpy.testing import (
    assert_allclose,
)
import os
import signal
import subprocess
import numpy as np
from .base import assert_allclose_with_logging
from pathlib import Path
import time
import logging
from .datafiles import NAMD_TOPOL, NAMD_CONF, NAMD_TRAJ, NAMD_PARAMS, NAMD_PSF
from .datafiles import NAMD_VEL, NAMD_FORCE
from importlib import resources
from pathlib import Path

PDBVELFACTOR=20.45482706

build_architecture = "openmpi"  # "multi-core", "gpu", "mpi"
sol_env_file = "openmpi-sol.sh"
namd_executable = "openmi-namd3"

# if build_architecture == "openmpi":
#     sol_env_file = "openmpi-sol.sh"
#     namd_executable = "openmi-namd3"
# if build_architecture == "mpi":
#     sol_env_file = "mpi-sol.sh"
#     namd_executable = "namd3"
# elif build_architecture == "multi-core":
#     sol_env_file = "sol.sh"
#     namd_executable = "namd3"
# elif build_architecture == "gpu":
#     sol_env_file = "gpu-sol.sh"
#     namd_executable = "gpu-namd3"
# else:
#     print(build_architecture)
#     raise ValueError("Unsupported build architecture")

num_nodes = 2
num_tasks = 2
num_tasks_per_node = 1
num_cores_per_task = 4
num_threads_per_core = 1
num_threads_per_task = num_cores_per_task * num_threads_per_core - 1 # extra core for NAMD internal communication and resource management

# make sure to run `srun --ntasks={num_tasks} --ntasks-per-node={num_tasks_per_node} --cpus-per-task={num_cores_per_task}` before running the test on it`

_data_ref = resources.files("imdclient.data")

NAMD_SOL_ENV_FILE = (_data_ref / "namd" / "md" / sol_env_file).as_posix()
NAMD_SOL_ENV_CMD = "bash " + NAMD_SOL_ENV_FILE
NAMD_SOL_CMD = "mpiexec -n {num_tasks}"
NAMD_EXEC = (_data_ref / "namd" / "md" / namd_executable).as_posix() + f" ++ppn {num_cores_per_task}"

NAMD_IMD_TRAJ = "alanin_imd.trr"

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("sol_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3NAMDsol:
    @pytest.fixture()
    def true_u_traj(self):
        return mda.Universe(
            NAMD_TOPOL,
            NAMD_TRAJ,
        )

    @pytest.fixture()
    def true_u_vel(self):
        return mda.Universe(
            NAMD_TOPOL,
            NAMD_VEL,
        )
    
    @pytest.fixture()
    def true_u_force(self):
        return mda.Universe(
            NAMD_TOPOL,
            NAMD_FORCE,
        )
    
    @pytest.fixture()
    def imd_u(self, run_sim_and_wait, tmp_path):
        tmp_u = mda.Universe(
            NAMD_TOPOL,
            "imd://localhost:8888",
        )
        with mda.Writer(
            f"{tmp_path.as_posix()}/" + NAMD_IMD_TRAJ, tmp_u.atoms.n_atoms
        ) as w:
            for ts in tmp_u.trajectory:
                w.write(tmp_u.atoms)
        imd_u = mda.Universe(
            NAMD_TOPOL,
            f"{tmp_path.as_posix()}/" + NAMD_IMD_TRAJ,
        )
        # Give all MPI ranks a chance to release FD on traj
        time.sleep(10)
        yield imd_u

    @pytest.fixture()
    def command(self):
        return (
            f"cp {NAMD_PARAMS} {NAMD_PSF} {NAMD_TOPOL} . && {NAMD_SOL_ENV_CMD} && {NAMD_SOL_CMD} {NAMD_EXEC} {NAMD_CONF}"
        )

    @pytest.fixture()
    def run_sim_and_wait(self, tmp_path, command, match_string):
        old_cwd = Path.cwd()
        os.chdir(tmp_path)
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
        os.chdir(old_cwd)
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)

    @pytest.fixture()
    def match_string(self):
        return "INTERACTIVE MD AWAITING CONNECTION"

    @pytest.fixture()
    def first_frame(self):
        return 0

    @pytest.fixture()
    def universe(self):
        return mda.Universe(
            NAMD_TOPOL,
            NAMD_TRAJ,
        )

    def test_compare_imd_to_true_traj(self, true_u_traj, true_u_vel, true_u_force, imd_u, first_frame):

        for i in range(first_frame, len(true_u_traj.trajectory)):
            assert_allclose(
                true_u_traj.trajectory[i].time,
                imd_u.trajectory[i - first_frame].time,
                atol=1e-03,
            )
            assert_allclose(
                true_u_traj.trajectory[i].data["step"],
                imd_u.trajectory[i - first_frame].data["step"],
            )
            if true_u_traj.trajectory[i].dimensions is not None:
                assert_allclose_with_logging(
                    true_u_traj.trajectory[i].dimensions,
                    imd_u.trajectory[i - first_frame].dimensions,
                    atol=1e-03,
                )
            if true_u_traj.trajectory[i].has_positions:
                assert_allclose_with_logging(
                    true_u_traj.trajectory[i].positions,
                    imd_u.trajectory[i - first_frame].positions,
                    atol=1e-03,
                )

        for i in range(first_frame, len(true_u_vel.trajectory)):
            if true_u_vel.trajectory[i].has_positions:
                assert_allclose_with_logging(
                    true_u_vel.trajectory[i].positions*PDBVELFACTOR,
                    imd_u.trajectory[i - first_frame].positions,
                    atol=1e-03,
                )

        for i in range(first_frame, len(true_u_force.trajectory)):
            if true_u_force.trajectory[i].has_positions:
                assert_allclose_with_logging(
                    true_u_force.trajectory[i].positions,
                    imd_u.trajectory[i - first_frame].positions,
                    atol=1e-03,
                )
