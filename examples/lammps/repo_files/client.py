import MDAnalysis as mda
from MDAnalysis.coordinates.IMD import IMDReader
import logging
from imdclient.tests.datafiles import LAMMPS_TOPOL

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdreader.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

i = 0
u = mda.Universe(LAMMPS_TOPOL, "imd://localhost:8888")
for ts in u.trajectory:
    i += 1

logger.info(f"Parsed {i} frames")
