import struct
import logging
from enum import Enum, auto
from typing import Union
from dataclasses import dataclass
import numpy as np


IMDHEADERSIZE = 8
IMDENERGYPACKETLENGTH = 40
IMDBOXPACKETLENGTH = 36
IMDVERSIONS = {2, 3}
IMDAWAITGOTIME = 1

logger = logging.getLogger("imdclient.IMDClient")


class IMDHeaderType(Enum):
    IMD_DISCONNECT = 0
    IMD_ENERGIES = 1
    IMD_FCOORDS = 2
    IMD_GO = 3
    IMD_HANDSHAKE = 4
    IMD_KILL = 5
    IMD_MDCOMM = 6
    IMD_PAUSE = 7
    IMD_TRATE = 8
    IMD_IOERROR = 9
    # New in IMD v3
    IMD_SESSIONINFO = 10
    IMD_RESUME = 11
    IMD_TIME = 12
    IMD_BOX = 13
    IMD_VELOCITIES = 14
    IMD_FORCES = 15


def parse_energy_bytes(data, endianness):
    keys = [
        "step",
        "temperature",
        "total_energy",
        "potential_energy",
        "van_der_walls_energy",
        "coulomb_energy",
        "bonds_energy",
        "angles_energy",
        "dihedrals_energy",
        "improper_dihedrals_energy",
    ]
    values = struct.unpack(f"{endianness}ifffffffff", data)
    return dict(zip(keys, values))


def create_energy_bytes(
    step,
    temperature,
    total_energy,
    potential_energy,
    van_der_walls_energy,
    coulomb_energy,
    bonds_energy,
    angles_energy,
    dihedrals_energy,
    improper_dihedrals_energy,
    endianness,
):
    return struct.pack(
        f"{endianness}ifffffffff",
        step,
        temperature,
        total_energy,
        potential_energy,
        van_der_walls_energy,
        coulomb_energy,
        bonds_energy,
        angles_energy,
        dihedrals_energy,
        improper_dihedrals_energy,
    )


def parse_box_bytes(data, endianness):
    """Box is a 3x3 matrix of floats"""

    vals = struct.unpack(f"{endianness}fffffffff", data)
    return np.array(vals).reshape(3, 3)


class IMDHeader:
    """Convenience class to represent the header of an IMD packet"""

    def __init__(self, data):
        msg_type, length = struct.unpack("!ii", data)
        h_type = IMDHeaderType(msg_type)

        self.type = h_type
        self.length = length


class IMDTime:
    """Convenience class to represent the body of time packet"""

    def __init__(self, data, endianness):
        self.dt, self.time = struct.unpack(f"{endianness}ff", data)


@dataclass
class IMDSessionInfo:
    """Convenience class to represent the session information of an IMD connection

    '<' represents little endian and '>' represents big endian
    """

    version: int
    endianness: str
    # In IMDv2, we don't know if the server sends wrapped coordinates until
    # we receive the first packet.
    wrapped_coords: bool
    # In IMDv2, we don't know if the server sends energies until
    # we receive the first packet.
    energies: bool
    time: bool
    box: bool
    positions: bool
    velocities: bool
    forces: bool


def parse_imdv3_session_info(data, end):
    """Parses the session information packet of an IMD v3 connection"""
    logger.debug(f"parse_imdv3_session_info: {data}")
    time, energies, box, positions, velocties, forces, wrapped_coords = (
        struct.unpack(f"{end}BBBBBBB", data)
    )
    logger.debug(f"parse_imdv3_session_info2 : {data}")
    imdsinfo = IMDSessionInfo(
        version=3,
        endianness=end,
        time=(time == 1),
        energies=(energies == 1),
        box=(box == 1),
        positions=(positions == 1),
        velocities=(velocties == 1),
        forces=(forces == 1),
        wrapped_coords=(wrapped_coords == 1),
    )
    logger.debug(f"parse_imdv3_session_info3: {data}")
    return imdsinfo


def create_imdv3_session_info_bytes(imdsinfo: IMDSessionInfo):
    """Creates the session information packet of an IMD v3 connection"""
    body = 0
    body |= imdsinfo.time << 0
    body |= imdsinfo.energies << 1
    body |= imdsinfo.box << 2
    body |= imdsinfo.positions << 3
    body |= imdsinfo.velocities << 4
    body |= imdsinfo.forces << 5
    body |= imdsinfo.wrapped_coords << 6
    return struct.pack("!i", body)


def create_header_bytes(msg_type: IMDHeaderType, length: int):
    # NOTE: add error checking for invalid packet msg_type here

    type = msg_type.value
    return struct.pack("!ii", type, length)


def parse_header_bytes(data):
    msg_type, length = struct.unpack("!ii", data)
    type = IMDHeaderType(msg_type)
    # NOTE: add error checking for invalid packet msg_type here
    return IMDHeader(type, length)
