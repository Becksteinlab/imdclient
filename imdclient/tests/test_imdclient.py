"""Test for IMDClient functionality"""

from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_TRR,
    COORDINATES_H5MD,
)
import MDAnalysis as mda
import imdclient
from imdclient.IMDClient import imdframe_memsize, IMDClient
from imdclient.IMDProtocol import IMDHeaderType
from .utils import (
    get_free_port,
    create_default_imdsinfo_v2,
    create_default_imdsinfo_v3,
)
from .server import InThreadIMDServer
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseReference,
    BaseWriterTest,
    assert_timestep_almost_equal,
)
from MDAnalysisTests.coordinates.test_xdr import TRRReference
import numpy as np
import logging
import pytest
import time


logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)

IMDENERGYKEYS = [
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


class TestIMDClientV3:

    @pytest.fixture
    def port(self):
        return get_free_port()

    @pytest.fixture
    def universe(self):
        return mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)

    @pytest.fixture
    def imdsinfo(self):
        return create_default_imdsinfo_v3()

    @pytest.fixture
    def server_client(self, universe, imdsinfo, port):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", port, first_frame=False)
        client = IMDClient(
            f"localhost",
            port,
            universe.trajectory.n_atoms,
            buffer_size=imdframe_memsize(universe.trajectory.n_atoms, imdsinfo)
            * 2,
        )
        yield server, client
        client.stop()
        server.cleanup()

    def test_pause_resume_continue(self, server_client):
        server, client = server_client
        server.send_frames(0, 2)
        # Client's buffer is filled. client should send pause
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        # Empty buffer
        client.get_imdframe()
        # only the second call actually frees buffer memory
        client.get_imdframe()
        # client has free memory. should send resume
        server.expect_packet(IMDHeaderType.IMD_RESUME)
        server.send_frame(1)
        client.get_imdframe()

    def test_pause_resume_disconnect(self, server_client):
        """Client pauses because buffer is full, empties buffer and attempt to resume, but
        finds that simulation has already ended and raises EOF"""
        server, client = server_client
        server.send_frames(0, 2)
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        client.get_imdframe()
        client.get_imdframe()
        # client has free frame. should send resume
        server.expect_packet(IMDHeaderType.IMD_RESUME)
        # simulation is over. client should raise EOF
        server.disconnect()
        with pytest.raises(EOFError):
            client.get_imdframe()

    def test_pause_resume_no_disconnect(self, server_client):
        """Client pauses because buffer is full, empties buffer and attempt to resume, but
        finds that simulation has already ended (but has not yet disconnected) and raises EOF
        """
        server, client = server_client
        server.send_frames(0, 2)
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        client.get_imdframe()
        client.get_imdframe()
        # client has free frame. should send resume
        server.expect_packet(IMDHeaderType.IMD_RESUME)
        # simulation is over. client should raise EOF
        with pytest.raises(EOFError):
            client.get_imdframe()
        # server should receive disconnect from client (though it doesn't have to do anything)
        server.expect_packet(IMDHeaderType.IMD_DISCONNECT)


"""
class TestIMDClientV2:

    @pytest.fixture
    def port(self):
        return get_free_port()

    @pytest.fixture
    def traj(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def ref(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def server(self, traj):
        server = DummyIMDServer(traj, 2)
        return server

    @pytest.fixture(params=[">", "<"])
    def setup_test_endianness_traj_unchanged(self, request, server, port):
        server.port = port
        server.imdsessioninfo.endianness = request.param
        server.start()
        server.wait_for_event(IMDServerEventType.LISTENING)
        return server, port

    def test_endianness_traj_unchanged(
        self, setup_test_endianness_traj_unchanged, ref
    ):
        _, port = setup_test_endianness_traj_unchanged

        reader = imdreader.IMDREADER.IMDReader(
            f"localhost:{port}",
            convert_units=False,
            n_atoms=ref.trajectory.n_atoms,
        )

        i = 0
        # Can't call assert in loop- this prevents reader's __exit__ from being called
        # if assert fails. Instead copy timesteps and then assert them
        timesteps = []

        for ts in reader:
            logger.debug(
                f"test_imdreader: positions for frame {i}: {ts.positions}"
            )
            timesteps.append(ts.copy())
            i += 1

        assert i == len(ref)

        for j in range(len(ref)):
            np.testing.assert_allclose(timesteps[j].positions, ref[j].positions)
            offset = 0
            for energy_key in IMDENERGYKEYS:
                assert timesteps[j].data[energy_key] == j + offset
                offset += 1

    @pytest.fixture
    def setup_test_pause_traj_unchanged(self, server, port):
        server.port = port
        server.loop_behavior = ExpectPauseLoopV2Behavior()
        server.start()
        server.wait_for_event(IMDServerEventType.LISTENING)
        return server, port

    def test_pause_traj_unchanged(self, setup_test_pause_traj_unchanged, ref):
        server, port = setup_test_pause_traj_unchanged

        # Give the buffer only 1 IMDFrame of memory
        # We expect the producer thread to have to
        # pause every frame (except the first)
        reader = imdreader.IMDREADER.IMDReader(
            f"localhost:{port}",
            convert_units=False,
            n_atoms=ref.trajectory.n_atoms,
            buffer_size=imdframe_memsize(
                ref.trajectory.n_atoms, server.imdsessioninfo
            ),
        )

        i = 0
        timesteps = []

        for ts in reader:
            time.sleep(1)
            timesteps.append(ts.copy())
            i += 1

        assert i == len(ref)

        for j in range(len(ref)):
            np.testing.assert_allclose(timesteps[j].positions, ref[j].positions)
            offset = 0
            for energy_key in IMDENERGYKEYS:
                assert timesteps[j].data[energy_key] == j + offset
                offset += 1

    def test_no_connection(self):
        with pytest.raises(ConnectionRefusedError):
            imdreader.IMDREADER.IMDReader("localhost:12345", n_atoms=1)


class TestIMDReaderWithBlockingServerV2:

    @pytest.fixture
    def port(self):
        return get_free_port()

    @pytest.fixture
    def traj(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def ref(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def server(self, traj):
        server = TestIMDServer(traj)
        yield server
        server.cleanup()

    @pytest.mark.parametrize("endianness", [">", "<"])
    def test_change_endianness_traj_unchanged(self, ref, server, endianness):
        imdsinfo = create_default_imdsinfo_v2()
        imdsinfo.endianness = endianness

        host = "localhost"
        port = get_free_port()

        # This also sends first frame to prevent blocking
        # in reader's init
        server.listen_accept_handshake_send_ts(host, port, imdsinfo)
        reader = imdreader.IMDREADER.IMDReader(
            f"localhost:{port}",
            n_atoms=ref.trajectory.n_atoms,
            convert_units=False,
        )

        i = 0
        timesteps = []
        for ts in reader:
            if i != 4:
                server.send_frame(i + 1, endianness=endianness)
            if i == 4:
                server.disconnect()
            timesteps.append(ts.copy())
            i += 1

        assert i == len(ref)
"""
