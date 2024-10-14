from imdclient.IMDProtocol import *
import socket
from imdclient.IMDProtocol import *
import logging


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
