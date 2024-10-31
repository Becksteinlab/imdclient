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
LAMMPS_IN_NST_1 = (
    _data_ref / "lammps" / "md" / "lammps_v3_nst_1.in"
).as_posix()
LAMMPS_IN_NST_8 = (
    _data_ref / "lammps" / "md" / "lammps_v3_nst_8.in"
).as_posix()


GROMACS_GRO = (_data_ref / "gromacs" / "md" / "gromacs_struct.gro").as_posix()
GROMACS_MDP_NST_1 = (
    _data_ref / "gromacs" / "md" / "gromacs_v3_nst1.mdp"
).as_posix()
GROMACS_MDP_NST_8 = (
    _data_ref / "gromacs" / "md" / "gromacs_v3_nst8.mdp"
).as_posix()
GROMACS_TOP = (_data_ref / "gromacs" / "md" / "gromacs_v3.top").as_posix()

NAMD_TOPOL = (_data_ref / "namd" / "md" / "alanin.pdb").as_posix()
NAMD_CONF_NST_1 = (_data_ref / "namd" / "md" / "namd_v3_nst_1.namd").as_posix()
NAMD_CONF_NST_8 = (_data_ref / "namd" / "md" / "namd_v3_nst_8.namd").as_posix()
NAMD_PARAMS = (_data_ref / "namd" / "md" / "alanin.params").as_posix()
NAMD_PSF = (_data_ref / "namd" / "md" / "alanin.psf").as_posix()

del resources
