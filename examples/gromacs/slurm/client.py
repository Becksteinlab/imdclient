# NOTE: chamge
from imdclient.IMDREADER import IMDReader
import MDAnalysis as mda
import logging

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdreader.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

i = 0
u = mda.Universe("md.tpr", "localhost:8888")
for ts in u.trajectory:
    i += 1

logger.info(f"Parsed {i} frames")
