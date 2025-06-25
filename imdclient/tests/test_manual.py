from .minimalreader import MinimalReader
import MDAnalysis as mda
from numpy.testing import assert_allclose
import numpy as np
from pathlib import Path
import time
import argparse
import logging

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
    n_atoms = mda.Universe(
        topol_path,
        atom_style="id type x y z",
        convert_units=False,
    ).atoms.n_atoms
    tmp_u = MinimalReader(
        f"imd://localhost:8888", n_atoms=n_atoms, process_stream=True
    )
    return tmp_u


def test_compare_imd_to_true_traj_vel(imd_u, true_u_vel, first_frame):
    for i in range(first_frame, len(true_u_vel.trajectory)):
        # Manually convert unit
        assert_allclose_with_logging(
            true_u_vel.trajectory[i].positions * 20.45482706,
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


def test_compare_imd_to_true_traj(
    imd_u, true_u, first_frame, vel, force, dt, step
):
    for i in range(first_frame, len(true_u.trajectory)):
        assert_allclose(
            true_u.trajectory[i].time,
            imd_u.trajectory[i - first_frame].time,
            atol=1e-03,
        )
        # Issue #63
        # if dt:
        #     assert_allclose(
        #         true_u.trajectory[i].dt,
        #         imd_u.trajectory[i - first_frame].dt,
        #         atol=1e-03,
        #     )
        if step:
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
        if vel:
            assert_allclose_with_logging(
                true_u.trajectory[i].velocities,
                imd_u.trajectory[i - first_frame].velocities,
                atol=1e-03,
            )
        if force:
            assert_allclose_with_logging(
                true_u.trajectory[i].forces,
                imd_u.trajectory[i - first_frame].forces,
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

    print(
        "Writing IMD trajectory to temporary directory...\n===================="
    )
    imd_u = load_imd_universe(args.topol_path, args.tmp_path)

    print("Loading source of truth trajectory...\n====================")
    true_u = load_true_universe(args.topol_path, args.traj_path)

    try:
        print("Comparing trajectories...\n====================")
        vel_in_trr = args.vel_path is None
        force_in_trr = args.force_path is None
        dt_in_trr = not args.topol_path.endswith(".data")
        # True when not using DCDReader
        step_in_trr = not args.traj_path.endswith(".coor")

        test_compare_imd_to_true_traj(
            imd_u,
            true_u,
            args.first_frame,
            vel_in_trr,
            force_in_trr,
            dt_in_trr,
            step_in_trr,
        )

        if args.vel_path is not None:
            print(
                "Loading source of truth velocity trajectory...\n===================="
            )
            true_vel = load_true_universe(args.topol_path, args.vel_path)
            print("Comparing velocities...")
            test_compare_imd_to_true_traj_vel(imd_u, true_vel, args.first_frame)

        if args.force_path is not None:
            print(
                "Loading source of truth force trajectory...\n===================="
            )
            true_force = load_true_universe(args.topol_path, args.force_path)
            print("Comparing forces...")
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
