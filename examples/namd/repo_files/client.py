import MDAnalysis as mda
from imdclient.tests.minimalreader import minimalreader
import logging
from imdclient.tests.datafiles import NAMD_TOPOL

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdclient.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

i = 0
u = minimalreader("imd://localhost:8888")

while True:
    try:
        u._read_next_frame()
        i += 1
    except EOFError:
        break

logger.info(f"Parsed {i} frames")
