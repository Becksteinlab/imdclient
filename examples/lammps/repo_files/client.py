import logging

import MDAnalysis as mda

from imdclient.tests.minimalReader import minimalReader
from imdclient.tests.datafiles import LAMMPS_TOPOL

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdclient.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

i = 0
u_mda = mda.Universe(
    LAMMPS_TOPOL, atom_style="id type x y z", convert_units=False
)
n_atoms = u_mda.atoms.n_atoms
u = minimalReader("imd://localhost:8888", n_atoms=n_atoms)

while True:
    try:
        u._read_next_frame()
        i += 1
    except EOFError:
        break

logger.info(f"Parsed {i} frames")
