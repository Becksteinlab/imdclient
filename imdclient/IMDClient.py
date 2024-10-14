"""

IMDClient
^^^^^^^^^

.. autoclass:: IMDClient
   :members:

.. autoclass:: IMDProducerV3
   :members:
   :inherited-members:

.. autoclass:: IMDFrameBuffer
   :members:
   
"""

import socket
import threading
from .IMDProtocol import *
from .utils import read_into_buf, sock_contains_data, timeit
import logging
import queue
import time
import numpy as np
from typing import Union, Dict
import signal

logger = logging.getLogger(__name__)


class IMDClient:
    def __init__(
        self,
        host,
        port,
        n_atoms,
        socket_bufsize=None,
        multithreaded=True,
        **kwargs,
    ):
        """
        Parameters
        ----------
        host : str
            Hostname of the server
        port : int
            Port number of the server
        n_atoms : int
            Number of atoms in the simulation
        socket_bufsize : int, optional
            Size of the socket buffer in bytes. Default is to use the system default
        buffer_size : int, optional
            IMDFramebuffer will be filled with as many IMDFrames fit in `buffer_size` [``10MB``]
        **kwargs : optional
            Additional keyword arguments to pass to the IMDProducer and IMDFrameBuffer
        """
        self._stopped = False
        self._conn = self._connect_to_server(host, port, socket_bufsize)
        self._imdsinfo = self._await_IMD_handshake()
        self._multithreaded = multithreaded

        if self._multithreaded:
            self._buf = IMDFrameBuffer(
                self._imdsinfo,
                n_atoms,
                **kwargs,
            )
        else:
            self._buf = None
        if self._imdsinfo.version == 2:
            self._producer = IMDProducerV2(
                self._conn,
                self._buf,
                self._imdsinfo,
                n_atoms,
                multithreaded,
                **kwargs,
            )
        elif self._imdsinfo.version == 3:
            self._producer = IMDProducerV3(
                self._conn,
                self._buf,
                self._imdsinfo,
                n_atoms,
                multithreaded,
                **kwargs,
            )

        self._go()

        if self._multithreaded:
            signal.signal(signal.SIGINT, self.signal_handler)
            self._producer.start()

    def signal_handler(self, sig, frame):
        """Catch SIGINT to allow clean shutdown on CTRL+C
        This also ensures that main thread execution doesn't get stuck
        waiting in buf.pop_full_imdframe()"""
        self.stop()

    def get_imdframe(self):
        """
        Raises
        ------
        EOFError
            If there are no more frames to read from the stream
        """
        if self._multithreaded:
            try:
                return self._buf.pop_full_imdframe()
            except EOFError:
                # in this case, consumer is already finished
                # and doesn't need to be notified
                self._disconnect()
                self._stopped = True
                raise EOFError
        else:
            try:
                return self._producer._get_imdframe()
            except EOFError:
                self._disconnect()
                raise EOFError

    def get_imdsessioninfo(self):
        return self._imdsinfo

    def stop(self):
        if self._multithreaded:
            if not self._stopped:
                self._buf.notify_consumer_finished()
                self._disconnect()
                self._stopped = True
        else:
            self._disconnect()

    def _connect_to_server(self, host, port, socket_bufsize):
        """
        Establish connection with the server, failing out instantly if server is not running
        """
        conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if socket_bufsize is not None:
            # Default (linux):
            # /proc/sys/net/core/rmem_default
            # Max (linux):
            # /proc/sys/net/core/rmem_max
            conn.setsockopt(
                socket.SOL_SOCKET, socket.SO_RCVBUF, socket_bufsize
            )
        try:
            logger.debug(f"IMDClient: Connecting to {host}:{port}")
            conn.connect((host, port))
        except ConnectionRefusedError:
            raise ConnectionRefusedError(
                f"IMDClient: Connection to {host}:{port} refused"
            )
        return conn

    def _await_IMD_handshake(self) -> IMDSessionInfo:
        """
        Wait for the server to send a handshake packet, then determine
        IMD session information.
        """
        end = ">"
        ver = None

        self._conn.settimeout(5)

        h_buf = bytearray(IMDHEADERSIZE)
        try:
            read_into_buf(self._conn, h_buf)
        except (ConnectionError, TimeoutError, Exception) as e:
            logger.debug("IMDClient: No handshake packet received: %s", e)
            raise ConnectionError("IMDClient: No handshake packet received")

        header = IMDHeader(h_buf)

        if header.type != IMDHeaderType.IMD_HANDSHAKE:
            raise ValueError(
                f"Expected header type `IMD_HANDSHAKE`, got {header.type}"
            )

        if header.length not in IMDVERSIONS:
            # Try swapping endianness
            swapped = struct.unpack("<i", struct.pack(">i", header.length))[0]
            if swapped not in IMDVERSIONS:
                err_version = min(swapped, header.length)
                # NOTE: Add test for this
                raise ValueError(
                    f"Incompatible IMD version. Expected version in {IMDVERSIONS}, got {err_version}"
                )
            else:
                end = "<"
                ver = swapped
        else:
            ver = header.length

        sinfo = None
        if ver == 2:
            # IMD v2 does not send a configuration handshake body packet
            sinfo = IMDSessionInfo(
                version=ver,
                endianness=end,
                wrapped_coords=False,
                time=False,
                energies=True,
                box=False,
                positions=True,
                velocities=False,
                forces=False,
            )

        elif ver == 3:
            read_into_buf(self._conn, h_buf)
            header = IMDHeader(h_buf)
            if header.type != IMDHeaderType.IMD_SESSIONINFO:
                raise ValueError(
                    f"Expected header type `IMD_SESSIONINFO`, got {header.type}"
                )
            if header.length != 7:
                raise ValueError(
                    f"Expected header length 7, got {header.length}"
                )
            data = bytearray(7)
            read_into_buf(self._conn, data)
            sinfo = parse_imdv3_session_info(data, end)
            logger.debug(f"IMDClient: Received IMDv3 session info: {sinfo}")

        return sinfo

    def _go(self):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        go = create_header_bytes(IMDHeaderType.IMD_GO, 0)
        self._conn.sendall(go)
        logger.debug("IMDClient: Sent go packet to server")

    def _disconnect(self):
        # MUST disconnect before stopping execution
        # if simulation already ended, this method will do nothing
        try:
            disconnect = create_header_bytes(IMDHeaderType.IMD_DISCONNECT, 0)
            self._conn.sendall(disconnect)
            logger.debug("IMDClient: Disconnected from server")
        except (ConnectionResetError, BrokenPipeError, Exception) as e:
            logger.debug(
                f"IMDProducer: Attempted to disconnect but server already terminated the connection: %s",
                e,
            )
        finally:
            self._conn.close()


class BaseIMDProducer(threading.Thread):

    def __init__(
        self,
        conn,
        buffer,
        sinfo,
        n_atoms,
        multithreaded=True,
        timeout=5,
        **kwargs,
    ):
        """
        Parameters
        ----------
        conn : socket.socket
            Connection object to the server
        buffer : IMDFrameBuffer
            Buffer object to hold IMD frames. If `multithreaded` is False, this
            argument is ignored
        sinfo : IMDSessionInfo
            Information about the IMD session
        n_atoms : int
            Number of atoms in the simulation
        multithreaded : bool, optional
            If True, socket interaction will occur in a separate thread &
            frames will be buffered. Single-threaded, blocking IMDClient
            should only be used in testing [[``True``]]

        """
        super(BaseIMDProducer, self).__init__(daemon=True)
        self._conn = conn
        self._imdsinfo = sinfo
        self._paused = False

        # Timeout for first frame should be longer
        # than rest of frames
        self._timeout = timeout
        self._conn.settimeout(self._timeout)

        self._buf = buffer

        self._frame = 0
        self._parse_frame_time = 0

        self._n_atoms = n_atoms

        if not multithreaded:
            self._imdf = IMDFrame(n_atoms, sinfo)
        else:
            # In this case, buffer performs frame allocation
            self._imdf = None

        self._header = bytearray(IMDHEADERSIZE)

    def _parse_imdframe(self):
        raise NotImplementedError

    def _pause(self):
        raise NotImplementedError

    def _unpause(self):
        raise NotImplementedError

    def _get_imdframe(self):
        """Parses an IMD frame from the stream. Blocks until a frame is available.

        Raises
        ------
        EOFError
            If there are no more frames to read from the stream
        RuntimeError
            If an unexpected error occurs
        """
        try:
            self._parse_imdframe()
        except EOFError as e:
            raise EOFError
        except Exception as e:
            raise RuntimeError("An unexpected error occurred") from e

        return self._imdf

    def run(self):
        logger.debug("IMDProducer: Starting run loop")
        try:
            while True:
                with timeit() as t:
                    logger.debug("IMDProducer: Checking for pause")

                    if not self._paused:
                        logger.debug("IMDProducer: Not paused")
                        if self._buf.is_full():
                            self._pause()
                            self._paused = True

                    logger.debug("IMDProducer: Checking for unpause")

                    if self._paused:
                        logger.debug("IMDProducer: Paused")
                        # wait for socket to empty before unpausing
                        if not sock_contains_data(self._conn, 0):
                            logger.debug("IMDProducer: Checked sock for data")
                            self._buf.wait_for_space()
                            self._unpause()
                            self._paused = False

                    self._imdf = self._buf.pop_empty_imdframe()

                    logger.debug("IMDProducer: Got empty frame")

                    self._parse_imdframe()

                    self._buf.push_full_imdframe(self._imdf)

                self._frame += 1
                logger.debug(
                    "IMDProducer: Frame %d parsed in %d seconds",
                    self._frame,
                    t.elapsed,
                )
        except EOFError:
            # simulation ended in a way
            # that we expected
            # i.e. consumer stopped or read_into_buf didn't find
            # full token of data
            logger.debug("IMDProducer: Simulation ended normally, cleaning up")
        except Exception as e:
            logger.debug("IMDProducer: An unexpected error occurred: %s", e)
        finally:
            logger.debug("IMDProducer: Stopping run loop")
            # Tell consumer not to expect more frames to be added
            self._buf.notify_producer_finished()

    def _expect_header(self, expected_type, expected_value=None):

        header = self._get_header()

        if header.type != expected_type:
            raise RuntimeError(
                f"IMDProducer: Expected header type {expected_type}, got {header.type}"
            )
        # Sometimes we do not care what the value is
        if expected_value is not None and header.length != expected_value:
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


class IMDProducerV2(BaseIMDProducer):
    def __init__(self, conn, buffer, sinfo, n_atoms, multithreaded, **kwargs):
        super(IMDProducerV2, self).__init__(
            conn, buffer, sinfo, n_atoms, multithreaded, **kwargs
        )

        self._energies = bytearray(IMDENERGYPACKETLENGTH)
        self._positions = bytearray(12 * self._n_atoms)
        self._prev_energies = None

    def _parse_imdframe(self):
        # Not all IMDv2 implementations send energies
        # so we can't expect them

        # Even if they are sent, energies might not be sent every frame
        # cache the last energies received

        # Either receive energies + positions or just positions
        header = self._get_header()
        if header.type == IMDHeaderType.IMD_ENERGIES and header.length == 1:
            self._imdsinfo.energies = True
            self._read(self._energies)
            self._imdf.energies.update(
                parse_energy_bytes(self._energies, self._imdsinfo.endianness)
            )
            self._prev_energies = self._imdf.energies

            self._expect_header(
                IMDHeaderType.IMD_FCOORDS, expected_value=self._n_atoms
            )
            self._read(self._positions)
            np.copyto(
                self._imdf.positions,
                np.frombuffer(
                    self._positions, dtype=f"{self._imdsinfo.endianness}f"
                ).reshape((self._n_atoms, 3)),
            )
        elif (
            header.type == IMDHeaderType.IMD_FCOORDS
            and header.length == self._n_atoms
        ):
            # If we received positions but no energies
            # use the last energies received
            if self._prev_energies is not None:
                self._imdf.energies = self._prev_energies
            else:
                self._imdf.energies = None
            self._read(self._positions)
            np.copyto(
                self._imdf.positions,
                np.frombuffer(
                    self._positions, dtype=f"{self._imdsinfo.endianness}f"
                ).reshape((self._n_atoms, 3)),
            )
        else:
            raise RuntimeError("IMDProducer: Unexpected packet type or length")

    def _pause(self):
        self._conn.settimeout(0)
        logger.debug(
            "IMDProducer: Pausing simulation because buffer is almost full"
        )
        pause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(pause)
        except ConnectionResetError as e:
            # Simulation has already ended by the time we paused
            raise EOFError
        # Edge case: pause occured in the time between server sends its last frame
        # and closing socket
        # Simulation is not actually paused but is over, but we still want to read remaining data
        # from the socket

    def _unpause(self):
        self._conn.settimeout(self._timeout)
        logger.debug("IMDProducer: Unpausing simulation, buffer has space")
        unpause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(unpause)
        except ConnectionResetError as e:
            # Edge case: pause occured in the time between server sends its last frame
            # and closing socket
            # Simulation was never actually paused in this case and is now over
            raise EOFError
        # Edge case: pause & unpause occured in the time between server sends its last frame and closing socket
        # in this case, the simulation isn't actually unpaused but over


class IMDProducerV3(BaseIMDProducer):
    def __init__(
        self,
        conn,
        buffer,
        sinfo,
        n_atoms,
        multithreaded,
        **kwargs,
    ):
        super(IMDProducerV3, self).__init__(
            conn,
            buffer,
            sinfo,
            n_atoms,
            multithreaded,
            **kwargs,
        )
        # The body of an x/v/f packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        xvf_bytes = 12 * n_atoms
        if self._imdsinfo.energies:
            self._energies = bytearray(IMDENERGYPACKETLENGTH)
        if self._imdsinfo.time:
            self._time = bytearray(IMDTIMEPACKETLENGTH)
        if self._imdsinfo.box:
            self._box = bytearray(IMDBOXPACKETLENGTH)
        if self._imdsinfo.positions:
            self._positions = bytearray(xvf_bytes)
        if self._imdsinfo.velocities:
            self._velocities = bytearray(xvf_bytes)
        if self._imdsinfo.forces:
            self._forces = bytearray(xvf_bytes)

    def _pause(self):
        self._conn.settimeout(0)
        logger.debug(
            "IMDProducer: Pausing simulation because buffer is almost full"
        )
        pause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(pause)
        except ConnectionResetError as e:
            # Simulation has already ended by the time we paused
            raise EOFError
        # Edge case: pause occured in the time between server sends its last frame
        # and closing socket
        # Simulation is not actually paused but is over, but we still want to read remaining data
        # from the socket

    def _unpause(self):
        self._conn.settimeout(self._timeout)
        logger.debug("IMDProducer: Unpausing simulation, buffer has space")
        unpause = create_header_bytes(IMDHeaderType.IMD_RESUME, 0)
        try:
            self._conn.sendall(unpause)
        except ConnectionResetError as e:
            # Edge case: pause occured in the time between server sends its last frame
            # and closing socket
            # Simulation was never actually paused in this case and is now over
            raise EOFError
        # Edge case: pause & unpause occured in the time between server sends its last frame and closing socket
        # in this case, the simulation isn't actually unpaused but over

    def _parse_imdframe(self):
        if self._imdsinfo.time:
            self._expect_header(IMDHeaderType.IMD_TIME, expected_value=1)
            # use header buf to hold time data sicnce it is also
            # 8 bytes
            self._read(self._time)
            t = IMDTime(self._time, self._imdsinfo.endianness)
            self._imdf.dt = t.dt
            self._imdf.time = t.time
            self._imdf.step = t.step

            logger.debug(
                f"IMDProducer: Time: {self._imdf.time}, dt: {self._imdf.dt}, step: {self._imdf.step}"
            )
        if self._imdsinfo.energies:
            self._expect_header(IMDHeaderType.IMD_ENERGIES, expected_value=1)
            self._read(self._energies)
            self._imdf.energies.update(
                parse_energy_bytes(self._energies, self._imdsinfo.endianness)
            )
        if self._imdsinfo.box:
            self._expect_header(IMDHeaderType.IMD_BOX, expected_value=1)
            self._read(self._box)
            self._imdf.box = parse_box_bytes(
                self._box, self._imdsinfo.endianness
            )
        if self._imdsinfo.positions:
            self._expect_header(
                IMDHeaderType.IMD_FCOORDS, expected_value=self._n_atoms
            )
            self._read(self._positions)
            np.copyto(
                self._imdf.positions,
                np.frombuffer(
                    self._positions, dtype=f"{self._imdsinfo.endianness}f"
                ).reshape((self._n_atoms, 3)),
            )
        if self._imdsinfo.velocities:
            self._expect_header(
                IMDHeaderType.IMD_VELOCITIES, expected_value=self._n_atoms
            )
            self._read(self._velocities)
            np.copyto(
                self._imdf.velocities,
                np.frombuffer(
                    self._velocities, dtype=f"{self._imdsinfo.endianness}f"
                ).reshape((self._n_atoms, 3)),
            )
        if self._imdsinfo.forces:
            self._expect_header(
                IMDHeaderType.IMD_FORCES, expected_value=self._n_atoms
            )
            self._read(self._forces)
            np.copyto(
                self._imdf.forces,
                np.frombuffer(
                    self._forces, dtype=f"{self._imdsinfo.endianness}f"
                ).reshape((self._n_atoms, 3)),
            )

    def __del__(self):
        logger.debug("IMDProducer: I am being deleted")


class IMDFrameBuffer:
    """
    Acts as interface between producer (IMDProducer) and consumer (IMDClient) threads
    when IMDClient runs in multithreaded mode
    """

    def __init__(
        self,
        imdsinfo,
        n_atoms,
        buffer_size=(10 * 1024**2),
        pause_empty_proportion=0.25,
        unpause_empty_proportion=0.5,
        **kwargs,
    ):
        """
        Parameters
        ----------
        imdsinfo : IMDSessionInfo
            Information about the IMD session
        n_atoms : int
            Number of atoms in the simulation
        buffer_size : int, optional
            Size of the buffer in bytes [``10MB``]
        pause_empty_proportion : float, optional
            Lower threshold proportion of the buffer's IMDFrames that are empty
            before the simulation is paused [``0.25``]
        unpause_empty_proportion : float, optional
            Proportion of the buffer's IMDFrames that must be empty
            before the simulation is unpaused [``0.5``]
        """

        # Syncing reader and producer
        self._producer_finished = False
        self._consumer_finished = False

        self._prev_empty_imdf = None

        self._empty_q = queue.Queue()
        self._full_q = queue.Queue()
        self._empty_imdf_avail = threading.Condition(threading.Lock())
        self._full_imdf_avail = threading.Condition(threading.Lock())

        if pause_empty_proportion < 0 or pause_empty_proportion > 1:
            raise ValueError("pause_empty_proportion must be between 0 and 1")
        self._pause_empty_proportion = pause_empty_proportion
        if unpause_empty_proportion < 0 or unpause_empty_proportion > 1:
            raise ValueError(
                "unpause_empty_proportion must be between 0 and 1"
            )
        self._unpause_empty_proportion = unpause_empty_proportion

        if buffer_size <= 0:
            raise ValueError("Buffer size must be positive")
        # Allocate IMDFrames with all of xvf present in imdsinfo
        # even if they aren't sent every frame. Can be optimized if needed
        imdf_memsize = imdframe_memsize(n_atoms, imdsinfo)
        self._total_imdf = buffer_size // imdf_memsize
        logger.debug(
            f"IMDFrameBuffer: Total IMDFrames allocated: {self._total_imdf}"
        )
        if self._total_imdf == 0:
            raise ValueError(
                "Buffer size is too small to hold a single IMDFrame"
            )
        for i in range(self._total_imdf):
            self._empty_q.put(IMDFrame(n_atoms, imdsinfo))

        # Timing for analysis
        self._t1 = None
        self._t2 = None

        self._frame = 0

    def is_full(self):
        logger.debug("IMDFrameBuffer: Checking if full")
        if (
            self._empty_q.qsize() / self._total_imdf
            <= self._pause_empty_proportion
        ):
            return True

        logger.debug(
            f"IMDFrameBuffer: Full frames: {self._full_q.qsize()}/{self._total_imdf}"
        )
        return False

    def wait_for_space(self):
        logger.debug("IMDProducer: Waiting for space in buffer")

        # Before acquiring the lock, check if we can return immediately
        if (
            self._empty_q.qsize() / self._total_imdf
            >= self._unpause_empty_proportion
        ) and not self._consumer_finished:
            return
        try:
            with self._empty_imdf_avail:
                while (
                    self._empty_q.qsize() / self._total_imdf
                    < self._unpause_empty_proportion
                ) and not self._consumer_finished:
                    logger.debug("IMDProducer: Waiting...")
                    logger.debug(
                        f"IMDProducer: consumer_finished: {self._consumer_finished}"
                    )
                    self._empty_imdf_avail.wait()

            logger.debug(
                "IMDProducer: Got space in buffer or consumer finished"
            )

            if self._consumer_finished:
                logger.debug("IMDProducer: Noticing consumer finished")
                raise EOFError
        except Exception as e:
            logger.debug(
                f"IMDProducer: Error waiting for space in buffer: {e}"
            )

    def pop_empty_imdframe(self):
        logger.debug("IMDProducer: Getting empty frame")

        # If there are empty frames available, don't wait
        if self._empty_q.qsize() and not self._consumer_finished:
            return self._empty_q.get()

        # If there are no empty frames available, wait for one
        with self._empty_imdf_avail:
            while self._empty_q.qsize() == 0 and not self._consumer_finished:
                self._empty_imdf_avail.wait()

        if self._consumer_finished:
            logger.debug("IMDProducer: Noticing consumer finished")
            raise EOFError

        return self._empty_q.get()

    def push_full_imdframe(self, imdf):
        # If the full q is empty, the consumer is waiting
        # and needs to be awakened
        self._full_q.put(imdf)
        with self._full_imdf_avail:
            self._full_imdf_avail.notify()

    def pop_full_imdframe(self):
        """Put empty_ts in the empty_q and get the next full timestep"""
        # Start timer- one frame of analysis is starting (including removal
        # from buffer)
        self._t1 = self._t2
        self._t2 = time.time()

        # Return the processed IMDFrame
        if self._prev_empty_imdf is not None:
            self._empty_q.put(self._prev_empty_imdf)
            with self._empty_imdf_avail:
                self._empty_imdf_avail.notify()

        # Get the next IMDFrame
        imdf = None
        if self._full_q.qsize() and not self._consumer_finished:
            imdf = self._full_q.get()
        else:
            with self._full_imdf_avail:
                while (
                    self._full_q.qsize() == 0 and not self._producer_finished
                ):
                    self._full_imdf_avail.wait()

            if self._producer_finished and self._full_q.qsize() == 0:
                logger.debug("IMDFrameBuffer(Consumer): Producer finished")
                raise EOFError

            imdf = self._full_q.get()

        self._prev_empty_imdf = imdf

        if self._t1 is not None:
            logger.debug(
                f"IMDFrameBuffer(Consumer): Frame #{self._frame} analyzed in {self._t2 - self._t1} seconds"
            )
        self._frame += 1

        return imdf

    def notify_producer_finished(self):
        self._producer_finished = True
        with self._full_imdf_avail:
            self._full_imdf_avail.notify()

    def notify_consumer_finished(self):
        self._consumer_finished = True
        with self._empty_imdf_avail:
            # noop if producer isn't waiting
            self._empty_imdf_avail.notify()


class IMDFrame:

    def __init__(self, n_atoms, imdsinfo):
        if imdsinfo.time:
            self.time = 0.0
            self.dt = 0.0
            self.step = 0.0
        else:
            self.time = None
            self.dt = None
            self.step = None
        if imdsinfo.energies:
            self.energies = {
                "step": 0,
                "temperature": 0.0,
                "total_energy": 0.0,
                "potential_energy": 0.0,
                "van_der_walls_energy": 0.0,
                "coulomb_energy": 0.0,
                "bonds_energy": 0.0,
                "angles_energy": 0.0,
                "dihedrals_energy": 0.0,
                "improper_dihedrals_energy": 0.0,
            }
        else:
            self.energies = None
        if imdsinfo.box:
            self.box = np.empty((3, 3), dtype=np.float32)
        else:
            self.box = None
        if imdsinfo.positions:
            self.positions = np.empty((n_atoms, 3), dtype=np.float32)
        else:
            self.positions = None
        if imdsinfo.velocities:
            self.velocities = np.empty((n_atoms, 3), dtype=np.float32)
        else:
            self.velocities = None
        if imdsinfo.forces:
            self.forces = np.empty((n_atoms, 3), dtype=np.float32)
        else:
            self.forces = None


def imdframe_memsize(n_atoms, imdsinfo) -> int:
    """
    Calculate the memory size of an IMDFrame in bytes
    """
    memsize = 0
    if imdsinfo.time:
        memsize += 8 * 3
    if imdsinfo.energies:
        memsize += 4 * 10
    if imdsinfo.box:
        memsize += 4 * 9
    if imdsinfo.positions:
        memsize += 4 * 3 * n_atoms
    if imdsinfo.velocities:
        memsize += 4 * 3 * n_atoms
    if imdsinfo.forces:
        memsize += 4 * 3 * n_atoms

    return memsize
