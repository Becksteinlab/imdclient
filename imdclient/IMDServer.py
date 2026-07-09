import socket
from socket import SOL_SOCKET, SO_REUSEADDR
import threading
from .IMDProtocol import *
from .utils import read_into_buf, sock_contains_data, timeit
import logging
import queue
import time
import numpy as np
from typing import Union, Dict
import signal
import atexit
import sys

logger = logging.getLogger(__name__)


class IMDServer:

    def __init__(
        self,
        version,
        n_atoms,
        host="0.0.0.0",
        port=8888,
        time=True,
        box=True,
        positions=True,
        velocities=True,
        forces=True,
        listen_timeout=60,
        read_timeout=5,
        **kwargs,
    ):
        if version not in IMDVERSIONS:
            raise ValueError(
                f"IMDServer: Incompatible IMD version. Expected version in {IMDVERSIONS}, got {version}"
            )
        if version == 2 and (time or box or velocities or forces):
            # TODO: and energies?
            raise ValueError("IMDServer: IMDv2 only supports positions")
        if (version == 2 and not (positions)) or (
            version == 3
            and not (time or box or positions or velocities or forces)
        ):
            raise ValueError(
                f"IMDServer: No IMD flags turned on for IMDv{version}"
            )

        end = ">" if sys.byteorder == "big" else "<"

        self.sinfo = IMDSessionInfo(
            version=version,
            endianness=end,
            wrapped_coords=False,
            time=time,
            energies=False,
            box=box,
            positions=positions,
            velocities=velocities,
            forces=forces,
        )
        self.n_atoms = n_atoms

        # shutdown safety
        signal.signal(signal.SIGINT, self.signal_handler)
        signal.signal(signal.SIGTERM, self.signal_handler)
        try:
            import IPython
        except ImportError:
            has_ipython = False
        else:
            has_ipython = True

        if has_ipython:
            try:
                from IPython import get_ipython

                if get_ipython() is not None:
                    kernel = get_ipython().kernel
                    kernel.pre_handler_hook = lambda: None
                    kernel.post_handler_hook = lambda: None
                    logger.debug("Running in Jupyter")
            except NameError:
                logger.debug("Running in non-jupyter IPython environment")

        atexit.register(self.stop)

        self._alloc_send_buf()
        self._header = bytearray(IMDHEADERSIZE)
        self._array_dtype = np.dtype(f"{self.sinfo.endianness}f4")
        self._time_fmt = f"{self.sinfo.endianness}ddQ"

        try:
            self._listen_socket = socket.socket(
                socket.AF_INET, socket.SOCK_STREAM
            )
            self._listen_socket.setsockopt(SOL_SOCKET, SO_REUSEADDR, 1)
            self._listen_socket.bind((host, port))

            # for the case of port=0
            self.port = self._listen_socket.getsockname()[1]
            self._listen_socket.listen()

            if sock_contains_data(self._listen_socket, listen_timeout):
                self._conn, _ = self._listen_socket.accept()
            else:
                # TODO: abstract into loop
                raise TimeoutError(
                    f"IMDServer: No client within {listen_timeout} seconds"
                )

            if version == 2:
                self._send_handshakeV2()
            elif version == 3:
                self._send_handshakeV3()

            self._conn.settimeout(5)

            self._expect_header(IMDHeaderType.IMD_GO)

            self._conn.settimeout(read_timeout)
        except Exception as e:
            self.stop()
            raise e

    def write_frame(self, imdframe):

        if self.sinfo.time:
            if imdframe.time is None:
                raise ValueError(
                    f"IMDServer: Time enabled for session, but frame does not contain time information"
                )

            struct.pack_into(
                self._time_fmt,
                self._time_view,
                0,
                imdframe.dt,
                imdframe.time,
                imdframe.step,
            )

        if self.sinfo.box:
            if imdframe.box is None:
                raise ValueError(
                    f"IMDServer: Box enabled for session, but frame does not contain box information"
                )

            src = np.ascontiguousarray(imdframe.box, dtype=self._array_dtype)
            target = np.frombuffer(self._box_view, dtype=self._array_dtype)
            target = target.reshape(src.shape)
            np.copyto(target, src)

        if self.sinfo.positions:
            if imdframe.positions is None:
                raise ValueError(
                    f"IMDServer: Positions enabled for session, but frame does not contain positions information"
                )

            src = np.ascontiguousarray(
                imdframe.positions, dtype=self._array_dtype
            )
            target = np.frombuffer(self._pos_view, dtype=self._array_dtype)
            target = target.reshape(src.shape)
            np.copyto(target, src)

        if self.sinfo.velocities:
            if imdframe.velocities is None:
                raise ValueError(
                    f"IMDServer: Velocities enabled for session, but frame does not contain velocities information"
                )

            src = np.ascontiguousarray(
                imdframe.velocities, dtype=self._array_dtype
            )
            target = np.frombuffer(self._vel_view, dtype=self._array_dtype)
            target = target.reshape(src.shape)
            np.copyto(target, src)

        if self.sinfo.forces:
            if imdframe.forces is None:
                raise ValueError(
                    f"IMDServer: Forces enabled for session, but frame does not contain forces information"
                )

            src = np.ascontiguousarray(
                imdframe.forces, dtype=self._array_dtype
            )
            target = np.frombuffer(self._force_view, dtype=self._array_dtype)
            target = target.reshape(src.shape)
            np.copyto(target, src)

        self._conn.sendall(self._send_buf)

    def _alloc_send_buf(self):

        buf_len = 0
        xvf_size = 12 * self.n_atoms
        if self.sinfo.time:
            buf_len += IMDHEADERSIZE + IMDTIMEPACKETLENGTH
        if self.sinfo.box:
            buf_len += IMDHEADERSIZE + IMDBOXPACKETLENGTH
        if self.sinfo.positions:
            buf_len += IMDHEADERSIZE + xvf_size
        if self.sinfo.velocities:
            buf_len += IMDHEADERSIZE + xvf_size
        if self.sinfo.forces:
            buf_len += IMDHEADERSIZE + xvf_size

        self._send_buf = bytearray(buf_len)
        offset_bytes = 0

        if self.sinfo.time:
            self._send_buf[offset_bytes : offset_bytes + IMDHEADERSIZE] = (
                create_header_bytes(IMDHeaderType.IMD_TIME, 1)
            )
            offset_bytes += IMDHEADERSIZE
            self._time_view = memoryview(self._send_buf)[
                offset_bytes : offset_bytes + IMDTIMEPACKETLENGTH
            ]
            offset_bytes += IMDTIMEPACKETLENGTH

        if self.sinfo.box:
            self._send_buf[offset_bytes : offset_bytes + IMDHEADERSIZE] = (
                create_header_bytes(IMDHeaderType.IMD_BOX, 1)
            )
            offset_bytes += IMDHEADERSIZE
            self._box_view = memoryview(self._send_buf)[
                offset_bytes : offset_bytes + IMDBOXPACKETLENGTH
            ]
            offset_bytes += IMDBOXPACKETLENGTH

        if self.sinfo.positions:
            self._send_buf[offset_bytes : offset_bytes + IMDHEADERSIZE] = (
                create_header_bytes(IMDHeaderType.IMD_FCOORDS, self.n_atoms)
            )
            offset_bytes += IMDHEADERSIZE
            self._pos_view = memoryview(self._send_buf)[
                offset_bytes : offset_bytes + xvf_size
            ]
            offset_bytes += xvf_size

        if self.sinfo.velocities:
            self._send_buf[offset_bytes : offset_bytes + IMDHEADERSIZE] = (
                create_header_bytes(IMDHeaderType.IMD_VELOCITIES, self.n_atoms)
            )
            offset_bytes += IMDHEADERSIZE
            self._vel_view = memoryview(self._send_buf)[
                offset_bytes : offset_bytes + xvf_size
            ]
            offset_bytes += xvf_size

        if self.sinfo.forces:
            self._send_buf[offset_bytes : offset_bytes + IMDHEADERSIZE] = (
                create_header_bytes(IMDHeaderType.IMD_FORCES, self.n_atoms)
            )
            offset_bytes += IMDHEADERSIZE
            self._force_view = memoryview(self._send_buf)[
                offset_bytes : offset_bytes + xvf_size
            ]
            offset_bytes += xvf_size

    def _send_handshakeV2(self):
        header = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        header += struct.pack(f"{self.sinfo.endianness}i", 2)
        self._conn.sendall(header)

    def _send_handshakeV3(self):
        logger.debug(f"InThreadIMDServer: Sending handshake V3")
        packet = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        packet += struct.pack(f"{self.sinfo.endianness}i", 3)
        self._conn.sendall(packet)

        sinfo = struct.pack("!ii", IMDHeaderType.IMD_SESSIONINFO.value, 7)
        time = 1 if self.sinfo.time else 0
        energies = 1 if self.sinfo.energies else 0
        box = 1 if self.sinfo.box else 0
        positions = 1 if self.sinfo.positions else 0
        velocities = 1 if self.sinfo.velocities else 0
        forces = 1 if self.sinfo.forces else 0
        wrapped_coords = 0
        sinfo += struct.pack(
            f"{self.sinfo.endianness}BBBBBBB",
            time,
            energies,
            box,
            positions,
            wrapped_coords,
            velocities,
            forces,
        )
        logger.debug(f"IMDServer: Sending session info")
        self._conn.sendall(sinfo)

    def _expect_header(self, expected_type, expected_value=None):
        header = self._get_header()

        if header.type != expected_type:
            raise RuntimeError(
                f"IMDProducer: Expected header type {expected_type}, got {header.type}"
            )
        # Sometimes we do not care what the value is
        if expected_value is not None and header.length != expected_value:
            if expected_type in [
                IMDHeaderType.IMD_FCOORDS,
                IMDHeaderType.IMD_VELOCITIES,
                IMDHeaderType.IMD_FORCES,
            ]:
                raise RuntimeError(
                    f"IMDProducer: Expected n_atoms value {expected_value}, got {header.length}. "
                    + "Ensure you are using the correct topology file."
                )
            else:
                raise RuntimeError(
                    f"IMDProducer: Expected header value {expected_value}, got {header.length}"
                )

    def _get_header(self):
        self._read(self._header)
        return IMDHeader(self._header)

    def _read(self, buf):
        """Wraps `read_into_buf` call to give uniform error handling which indicates end of stream"""
        try:
            read_into_buf(self._conn, buf)
        except (ConnectionError, TimeoutError, BlockingIOError, Exception):
            # ConnectionError: Server is definitely done sending frames, socket is closed
            # TimeoutError: Server is *likely* done sending frames.
            # BlockingIOError: Occurs when timeout is 0 in place of a TimeoutError. Server is *likely* done sending frames
            # OSError: Occurs when main thread disconnects from the server and closes the socket, but producer thread attempts to read another frame
            # Exception: Something unexpected happened
            raise EOFError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()
        return False

    def stop(self):
        print("Stopping")
        try:
            self._conn.close()
        except:
            pass
        try:
            self._listen_socket.close()
        except:
            pass

    def signal_handler(self, *args, **kwargs):
        """Catch SIGINT to allow clean shutdown on CTRL+C."""
        logger.debug("Intercepted signal")
        self.stop()
        logger.debug("Shutdown success")
