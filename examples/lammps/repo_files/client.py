import MDAnalysis as mda
from imdclient.utils import parse_host_port
from imdclient.IMDClient import IMDClient
import logging
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
# Parse the host and port from the IMD producer server address
host, port = parse_host_port("imd://localhost:8888")

n_atoms = mda.Universe(
    LAMMPS_TOPOL, atom_style="id type x y z", convert_units=False
).atoms.n_atoms
# This starts the simulation
imdclient = IMDClient(host, port, n_atoms=n_atoms)

while True:
    try:
        imd_frame = imdclient.get_imdframe()
    except EOFError:
        break
    i += 1
    logger.debug(f"IMDClient: Received frame {i}")

logger.info(f"IMDClient: Parsed {i} frames")
