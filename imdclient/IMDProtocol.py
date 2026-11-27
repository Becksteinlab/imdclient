"""
IMDProtocol module
==================

Low-level IMDv2 / IMDv3 packet types, session metadata, and parsers.

.. autoclass:: IMDHeaderType
   :members:
   :undoc-members:

.. autoclass:: IMDHeader
   :members:

.. autoclass:: IMDTime
   :members:

.. autoclass:: IMDSessionInfo
   :members:

.. autofunction:: parse_energy_bytes
.. autofunction:: create_energy_bytes
.. autofunction:: parse_box_bytes
.. autofunction:: parse_imdv3_session_info
.. autofunction:: create_header_bytes
.. autofunction:: parse_header_bytes
"""

import struct
import logging
from enum import Enum, auto
from typing import Union
from dataclasses import dataclass
import numpy as np

IMDHEADERSIZE = 8
IMDENERGYPACKETLENGTH = 40
IMDBOXPACKETLENGTH = 36
IMDTIMEPACKETLENGTH = 24
IMDVERSIONS = {2, 3}
IMDAWAITGOTIME = 1

logger = logging.getLogger("imdclient.IMDClient")


class IMDHeaderType(Enum):
    """IMD packet header type codes (32-bit enum values when transmitted).

    Each member is the header type for one IMD packet. Authoritative packet
    layouts and semantics are in :doc:`protocol_v3`.
    """

    IMD_DISCONNECT = 0
    """Disconnect packet: receiver asks the engine to close the socket.

    See :ref:`disconnect`.
    """

    IMD_ENERGIES = 1
    """Energies packet: optional per-frame IMD energy block from the engine.

    See :ref:`energies`.
    """

    IMD_FCOORDS = 2
    """Coordinates packet: optional per-frame atomic positions from the engine.

    See :ref:`coordinates`.
    """

    IMD_GO = 3
    """Go packet: receiver signals it is ready; engine may start sending frames.

    See :ref:`go`.
    """

    IMD_HANDSHAKE = 4
    """Handshake packet: engine announces IMD version and native endianness.

    See :ref:`handshake`.
    """

    IMD_KILL = 5
    """Kill packet: receiver requests that the engine stop and exit.

    See :ref:`kill`.
    """

    IMD_MDCOMM = 6
    """MD communication packet: receiver sends forces to apply to IMD atoms.

    See :ref:`md-communication`.
    """

    IMD_PAUSE = 7
    """Pause packet: receiver requests that the engine pause until resume.

    See :ref:`pause`.
    """

    IMD_TRATE = 8
    """Transmission-rate packet: receiver sets the IMD frame spacing.

    See :ref:`transmission-rate`.
    """

    IMD_IOERROR = 9
    """IO error code: not exchanged on the wire; for internal error signaling.

    See :ref:`io-error`.
    """

    # New in IMD v3
    IMD_SESSIONINFO = 10
    """Session-info packet (IMDv3): which per-frame data the engine will send.

    See :ref:`session-info`.
    """

    IMD_RESUME = 11
    """Resume packet (IMDv3): receiver requests resume after pause.

    See :ref:`resume`.
    """

    IMD_TIME = 12
    """Time packet (IMDv3): optional per-frame ``dt``, time, and step.

    See :ref:`time`.
    """

    IMD_BOX = 13
    """Box packet (IMDv3): optional per-frame triclinic box vectors.

    See :ref:`box`.
    """

    IMD_VELOCITIES = 14
    """Velocities packet (IMDv3): optional per-frame atomic velocities.

    See :ref:`velocities`.
    """

    IMD_FORCES = 15
    """Forces packet (IMDv3): optional per-frame atomic forces from the engine.

    See :ref:`forces`.
    """

    IMD_WAIT = 16
    """Wait packet (IMDv3): receiver sets blocking vs non-blocking wait behavior.

    See :ref:`wait`.
    """


def parse_energy_bytes(data, endianness):
    """Parse an IMD energies body into a dict of named fields.

    Parameters
    ----------
    data : bytes or bytearray
        Energy body payload (``int32`` step + nine ``float32`` values).
    endianness : str
        ``'<'`` or ``'>'``; native endianness of the simulation engine
        (from the handshake).

    Returns
    -------
    dict
        Keys ``step``, ``temperature``, ``total_energy``,
        ``potential_energy``, ``van_der_walls_energy``, ``coulomb_energy``,
        ``bonds_energy``, ``angles_energy``, ``dihedrals_energy``,
        ``improper_dihedrals_energy``. Same keys as
        :attr:`imdclient.IMDClient.IMDFrame.energies`.
    """
    keys = [
        "step",
        "temperature",
        "total_energy",
        "potential_energy",
        "van_der_waals_energy",
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
    van_der_waals_energy,
    coulomb_energy,
    bonds_energy,
    angles_energy,
    dihedrals_energy,
    improper_dihedrals_energy,
    endianness,
):
    """Pack named energy fields into an IMD energies body.

    Parameters
    ----------
    step : int
        Integration step (``int32``).
    temperature : float
        Absolute temperature (``float32``).
    total_energy, potential_energy, van_der_walls_energy, coulomb_energy, bonds_energy, angles_energy, dihedrals_energy, improper_dihedrals_energy : float
        Remaining energy-block ``float32`` terms (same names as
        :func:`parse_energy_bytes`).
    endianness : str
        ``'<'`` or ``'>'``; engine-native byte order.

    Returns
    -------
    bytes
        Packed energies body suitable to send after an
        :attr:`~imdclient.IMDProtocol.IMDHeaderType.IMD_ENERGIES` header.
    """
    return struct.pack(
        f"{endianness}ifffffffff",
        step,
        temperature,
        total_energy,
        potential_energy,
        van_der_waals_energy,
        coulomb_energy,
        bonds_energy,
        angles_energy,
        dihedrals_energy,
        improper_dihedrals_energy,
    )


def parse_box_bytes(data, endianness):
    """Parse an IMD box body into triclinic box vectors.

    Parameters
    ----------
    data : bytes or bytearray
        Box body payload (nine ``float32`` values).
    endianness : str
        ``'<'`` or ``'>'``; native endianness of the simulation engine
        (from the handshake).

    Returns
    -------
    ndarray
        Shape ``(3, 3)`` array of box vectors (rows are cell vectors).
    """
    vals = struct.unpack(f"{endianness}fffffffff", data)
    return np.array(vals).reshape(3, 3)


class IMDHeader:
    """Parsed 8-byte IMD packet header (network / big endian).

    Parameters
    ----------
    data : bytes or bytearray
        Exactly eight header bytes.

    Attributes
    ----------
    type : IMDHeaderType
        Packet header type.
    length : int
        Header length / data slot (``int32``). Meaning depends on
        :attr:`type` (unused for some packets; atom count, count of
        blocks, or a flag for others).
    """

    def __init__(self, data):
        msg_type, length = struct.unpack("!ii", data)
        h_type = IMDHeaderType(msg_type)

        self.type = h_type
        self.length = length


class IMDTime:
    """Parsed IMDv3 time-packet body (``dt``, time, step).

    Parameters
    ----------
    data : bytes or bytearray
        Time body payload (``float64`` dt, ``float64`` time, ``int64`` step).
    endianness : str
        ``'<'`` or ``'>'``; engine-native byte order from the handshake.

    Attributes
    ----------
    dt : float
        Integrator time step in ps.
    time : float
        Current simulation time in ps.
    step : int
        Current integration step (64-bit on the wire). Prefer this over
        ``energies["step"]`` for IMDv3 analysis.
    """

    def __init__(self, data, endianness):
        self.dt, self.time, self.step = struct.unpack(f"{endianness}ddq", data)


@dataclass
class IMDSessionInfo:
    """Negotiated IMD session metadata from the handshake / session-info.

    For IMDv3 this is filled from the session-info body. For IMDv2, several
    flags are inferred or remain provisional until the first frame
    (energies, wrapped coordinates).

    Attributes
    ----------
    version : int
        IMD protocol version for the session (``2`` or ``3``).
    endianness : str
        ``'<'`` (little) or ``'>'`` (big); native endianness of the engine,
        taken from the handshake and used when unpacking body packets.
    wrapped_coords : bool
        If ``True``, coordinates are wrapped into the simulation box.
        For IMDv2 this may be unknown until the first frame.
    energies : bool
        If ``True``, each IMD frame includes an energies packet.
        For IMDv2 this may be provisional until the first frame.
    time : bool
        If ``True``, each IMD frame includes a time packet (IMDv3).
    box : bool
        If ``True``, each IMD frame includes a box packet (IMDv3).
    positions : bool
        If ``True``, each IMD frame includes coordinates.
    velocities : bool
        If ``True``, each IMD frame includes velocities (IMDv3).
    forces : bool
        If ``True``, each IMD frame includes forces (IMDv3).
    """

    version: int
    endianness: str
    wrapped_coords: bool
    energies: bool
    time: bool
    box: bool
    positions: bool
    velocities: bool
    forces: bool


def parse_imdv3_session_info(data, end):
    """Parse an IMDv3 session-info body into :class:`IMDSessionInfo`.

    Parameters
    ----------
    data : bytes or bytearray
        Session-info body (seven ``uint8`` flags: time, energies, box,
        positions, wrapped_coords, velocities, forces).
    end : str
        ``'<'`` or ``'>'``; engine endianness from the handshake. Stored on
        the returned :class:`IMDSessionInfo` and used for later body parsing.

    Returns
    -------
    IMDSessionInfo
        Session flags with ``version=3`` and nonzero body bytes mapped to
        ``True``.
    """
    logger.debug(f"parse_imdv3_session_info: {data}")
    time, energies, box, positions, wrapped_coords, velocties, forces = (
        struct.unpack(f"{end}BBBBBBB", data)
    )
    logger.debug(f"parse_imdv3_session_info2 : {data}")
    imdsinfo = IMDSessionInfo(
        version=3,
        endianness=end,
        time=(time != 0),
        box=(box != 0),
        positions=(positions != 0),
        wrapped_coords=(wrapped_coords != 0),
        velocities=(velocties != 0),
        forces=(forces != 0),
        energies=(energies != 0),
    )
    logger.debug(f"parse_imdv3_session_info3: {data}")
    return imdsinfo


def create_header_bytes(msg_type: IMDHeaderType, length: int):
    """Pack an 8-byte IMD header (network / big endian).

    Parameters
    ----------
    msg_type : IMDHeaderType
        Packet header type to encode.
    length : int
        Header length / data slot (``int32``); meaning depends on
        ``msg_type`` (unused for some packets, atom count or flag for
        others).

    Returns
    -------
    bytes
        Eight-byte header ``(!ii)`` ready to send on the socket.
    """
    # NOTE: add error checking for invalid packet msg_type here

    type = msg_type.value
    return struct.pack("!ii", type, length)


def parse_header_bytes(data):
    """Parse an 8-byte IMD header from network / big-endian bytes.

    Parameters
    ----------
    data : bytes or bytearray
        Exactly eight bytes: type and length as ``int32`` each.

    Returns
    -------
    IMDHeader
        Header with :attr:`~IMDHeader.type` as :class:`IMDHeaderType` and
        :attr:`~IMDHeader.length` as ``int``.
    """
    msg_type, length = struct.unpack("!ii", data)
    type = IMDHeaderType(msg_type)
    # NOTE: add error checking for invalid packet msg_type here
    return IMDHeader(type, length)
