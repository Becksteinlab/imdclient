import imdclient
from imdclient.IMD import IMDReader
import MDAnalysis as mda
from numpy.testing import assert_allclose
import numpy as np
from pathlib import Path
import time
import argparse
import logging
from base import assert_allclose_with_logging

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("manual_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)

"""
Tool for running IMDv3 integration tests via the command line.

To use, start the simulation, wait for it to be ready for an IMD connection,
and then run this command relative to the root of the cloned respository:

python imdclient/tests/test_manual.py \
    --topol_path <path/to/topology> \
    --traj_path <path/to/trajectory> \
    --first_frame <first frame to compare> \
    --tmp_path <temporary directory for IMD trajectory>

Where the topology is the same topology as the IMD system, the trajectory is the path where
the trajectory of the running simulation is being written, and the first frame is the first frame of the
trajectory which should be compared to IMD data read from the socket (0 for GROMACS and NAMD, 1 for LAMMPS)
"""


def load_true_universe(topol_path, traj_path):
    if topol_path.endswith(".data"):
        return mda.Universe(
            topol_path,
            traj_path,
            atom_style="id type x y z",
            convert_units=False,
        )
    return mda.Universe(topol_path, traj_path)


def load_imd_universe(topol_path, tmp_path):
    # Pass atom_style (ignored if not using LAMMPS topol)
    tmp_u = mda.Universe(
        topol_path,
        "imd://localhost:8888",
        atom_style="id type x y z",
    )
    tmp_traj_file = f"{tmp_path}/imd_test_traj.trr"
    with mda.Writer(tmp_traj_file, tmp_u.atoms.n_atoms) as w:
        for ts in tmp_u.trajectory:
            w.write(tmp_u.atoms)
    time.sleep(10)  # Give MPI ranks a chance to release FD
    return mda.Universe(topol_path, tmp_traj_file, atom_style="id type x y z")


def test_compare_imd_to_true_traj_vel(imd_u, true_u_vel, first_frame):
    for i in range(first_frame, len(true_u_vel.trajectory)):
        assert_allclose_with_logging(
            true_u_vel.trajectory[i].positions,
            imd_u.trajectory[i - first_frame].velocities,
            atol=1e-03,
        )


def test_compare_imd_to_true_traj_forces(imd_u, true_u_force, first_frame):
    for i in range(first_frame, len(true_u_force.trajectory)):
        assert_allclose_with_logging(
            true_u_force.trajectory[i].positions,
            imd_u.trajectory[i - first_frame].forces,
            atol=1e-03,
        )


def test_compare_imd_to_true_traj(imd_u, true_u, first_frame):
    for i in range(first_frame, len(true_u.trajectory)):
        assert_allclose(
            true_u.trajectory[i].time,
            imd_u.trajectory[i - first_frame].time,
            atol=1e-03,
        )
        assert_allclose(
            true_u.trajectory[i].dt,
            imd_u.trajectory[i - first_frame].dt,
            atol=1e-03,
        )
        assert_allclose(
            true_u.trajectory[i].data["step"],
            imd_u.trajectory[i - first_frame].data["step"],
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


def main():
    parser = argparse.ArgumentParser(description="IMDv3 Integration Test Tool")
    parser.add_argument(
        "--topol_path", required=True, help="Path to topology file"
    )
    parser.add_argument(
        "--traj_path", required=True, help="Path to trajectory file"
    )
    parser.add_argument(
        "--vel_path",
        required=False,
        help="Path to velocities trajectory file (NAMD only)",
    )
    parser.add_argument(
        "--force_path",
        required=False,
        help="Path to forces trajectory file (NAMD only)",
    )
    parser.add_argument(
        "--first_frame",
        type=int,
        required=True,
        help="First frame to compare (0 for GROMACS and NAMD, 1 for LAMMPS)",
    )
    parser.add_argument(
        "--tmp_path", default=".", help="Temporary directory for IMD traj"
    )

    args = parser.parse_args()

    print("Writing IMD trajectory to temporary directory...")
    imd_u = load_imd_universe(args.topol_path, args.tmp_path)

    print("Loading source of truth trajectory...")
    true_u = load_true_universe(args.topol_path, args.traj_path)

    try:
        print("Comparing trajectories...")
        test_compare_imd_to_true_traj(true_u, imd_u, args.first_frame)

        if args.vel_path is not None:
            print("Loading source of truth velocity trajectory...")
            true_vel = load_true_universe(args.topol_path, args.vel_path)
            print("Comparing velocities...")
            test_compare_imd_to_true_traj_vel(
                imd_u, true_vel, args.first_frame
            )

        if args.force_path is not None:
            logger.info("Loading source of truth force trajectory...")
            true_force = load_true_universe(args.topol_path, args.force_path)
            logger.info("Comparing forces...")
            test_compare_imd_to_true_traj_forces(
                imd_u, true_force, args.first_frame
            )

        print("All tests passed!")

    except AssertionError as e:
        logger.error("Comparison failed!")
        print(f"Test failed: {e}")
        raise e


if __name__ == "__main__":
    main()
