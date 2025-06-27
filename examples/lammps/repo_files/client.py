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

# Parse the host and port from the IMD producer server address
host, port = parse_host_port("imd://localhost:8888")

n_atoms = mda.Universe(
    LAMMPS_TOPOL, atom_style="id type x y z", convert_units=False
).atoms.n_atoms
# This starts the simulation
imdclient = IMDClient(host, port, n_atoms=n_atoms)

atom_index = 0

i = 0
while True:
    try:
        imd_frame = imdclient.get_imdframe()
    except EOFError:
        break
    else:
        i += 1
        # `imd_frame` is an IMDFrame object containing frame data
        # Below is an example of how to access the data
        # For example, print frame number, simulation time, and positions of atom at `atom_index`
        print(
            f"Frame {i}: time={imd_frame.time}, atom {atom_index} position={imd_frame.positions[atom_index]}"
        )
        # You can also access other attributes like energies, box dimensions, velocities and forces,
        logger.debug(f"IMDClient: Received frame {i}")

logger.info(f"IMDClient: Parsed {i} frames")
