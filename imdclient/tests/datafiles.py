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

LAMMPS_TOPOL = (_data_ref / "lammps" / "md" / "lammps_topol.data").as_posix()
LAMMPS_IN = (_data_ref / "lammps" / "md" / "lammps_v3.in").as_posix()
LAMMPS_TRAJ = (_data_ref / "lammps" / "md" / "lammps_trj.h5md").as_posix()
GROMACS_TRAJ = (
    _data_ref / "gromacs" / "md" / "gromacs_v3_nst1.trr"
).as_posix()
GROMACS_TOPOL = (
    _data_ref / "gromacs" / "md" / "gromacs_struct.gro"
).as_posix()
GROMACS_TPR = (_data_ref / "gromacs" / "md" / "gromacs_v3_nst1.tpr").as_posix()
NAMD_TOPOL = (_data_ref / "namd" / "md" / "alanin.pdb").as_posix()
NAMD_CONF = (_data_ref / "namd" / "md" / "namd_v3.namd").as_posix()
NAMD_TRAJ = (_data_ref / "namd" / "md" / "alanin.dcd").as_posix()
NAMD_PARAMS = (_data_ref / "namd" / "md" / "alanin.params").as_posix()
NAMD_PSF = (_data_ref / "namd" / "md" / "alanin.psf").as_posix()

del resources
