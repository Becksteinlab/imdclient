"""
Location of data files
======================

Use as ::

    from imdclient.tests.datafiles import *

"""

__all__ = ["LAMMPS_IN", "LAMMPS_TOPOL", "GROMACS_TRAJ", "GROMACS_MDP"]

from importlib import resources
from pathlib import Path

_data_ref = resources.files("imdclient.data")

LAMMPS_TOPOL = (_data_ref / "lammps_topol.data").as_posix()
LAMMPS_IN = (_data_ref / "lammps_v3.in").as_posix()
LAMMPS_TRAJ = (_data_ref / "lammps_trj.h5md").as_posix()
GROMACS_TRAJ = (_data_ref / "gromacs_trj.trr").as_posix()
GROMACS_TOPOL = (_data_ref / "gromacs_topol.tpr").as_posix()


del resources
