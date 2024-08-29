from imdclient.IMDProtocol import *
import socket
import time
from imdclient.IMDProtocol import *
from imdclient.IMDREADER import read_into_buf, sock_contains_data
import logging
import subprocess

logger = logging.getLogger("imdclient.IMDClient")


def create_default_imdsinfo_v2():
    return IMDSessionInfo(
        version=2,
        endianness="<",
        wrapped_coords=True,
        energies=True,
        box=False,
        positions=True,
        velocities=False,
        forces=False,
    )


def create_default_imdsinfo_v3():
    return IMDSessionInfo(
        version=3,
        endianness="<",
        time=True,
        energies=True,
        box=True,
        positions=True,
        velocities=True,
        forces=True,
        wrapped_coords=False,
    )


def get_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def run_sim_and_wait(command, match_string):
    p = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        text=True,
        bufsize=0,
    )
    t = time.time()
    for stdout_line in iter(p.stdout.readline, ""):
        logger.debug(f"stdout: {stdout_line}")
        if match_string in stdout_line:
            break
        if time.time() - t > 10:
            raise TimeoutError("Timeout waiting for match string")

    logger.debug("Match string found")
    return p
