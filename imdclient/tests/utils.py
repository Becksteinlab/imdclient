from imdclient.IMDProtocol import *
import socket
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

def parse_host_port(filename):
    if not filename.startswith("imd://"):
        raise ValueError("IMDClient: URL must be in the format 'imd://host:port'")
    
    # Check if the format is correct
    parts = filename.split("imd://")[1].split(":")
    if len(parts) == 2:
        host = parts[0] 
        try:
            port = int(parts[1])
            return (host, port)
        except ValueError:
            raise ValueError("IMDClient: Port must be an integer")
    else:
        raise ValueError("IMDClient: URL must be in the format 'imd://host:port'")
