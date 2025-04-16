from MDAnalysis.coordinates.IMD import IMDReader
import MDAnalysis as mda
import logging
from imdclient.tests.datafiles import NAMD_TOPOL
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os

# Set up the FFmpeg command for streaming
# ffmpeg_command = [
#     'ffmpeg',
#     '-y',  # Overwrite output file if it exists
#     '-f', 'rawvideo',  # Input format
#     '-vcodec', 'rawvideo',  # Input codec
#     '-pix_fmt', 'rgb24',  # Pixel format
#     '-s', '640x480',  # Frame size
#     '-r', '1',  # Frame rate
#     '-i', '-',  # Input from standard input
#     '-an',  # No audio
#     '-c:v', 'libx264',  # Output codec
#     '-f', 'flv',  # Output format (e.g., flv for RTMP streaming)
#     'rtmp://localhost/live/stream'  # Output URL (or file)
# ]

# # Start FFmpeg subprocess
# ffmpeg = subprocess.Popen(ffmpeg_command, stdin=subprocess.PIPE)

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdclient.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

i = 0
u = mda.Universe(NAMD_TOPOL, "localhost:8888")
for ts in u.trajectory:
    # print(ts.positions[1][:2])
    i += 1
    
    fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=100)
    ax.plot(ts.positions[1][0], ts.positions[1][1], 'ro')
    ax.set_title(f"Position of atom 1")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    fig.savefig("Live-frame.png")
    
    # Render plot and convert it to raw RGB data
    # fig.canvas.draw()
    # frame = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    # frame = frame.reshape((480, 640, 3))  # Make sure this matches the size in ffmpeg_command
    # plt.close(fig)

    # # Write frame to FFmpeg's stdin as raw data
    # ffmpeg.stdin.write(frame.tobytes())

# ffmpeg.stdin.close()
# ffmpeg.wait()
logger.info(f"Parsed {i} frames")
