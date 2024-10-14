"""Test for IMDClient functionality"""

from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_H5MD,
)
import MDAnalysis as mda
from imdclient.IMDClient import imdframe_memsize, IMDClient
from imdclient.IMDProtocol import IMDHeaderType
from .utils import (
    get_free_port,
    create_default_imdsinfo_v3,
)
from .server import InThreadIMDServer
from MDAnalysisTests.coordinates.base import (
    assert_allclose,
)
from MDAnalysisTests.coordinates.test_xdr import TRRReference
import logging
import pytest


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
    def server_client_two_frame_buf(self, universe, imdsinfo, port):
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

    @pytest.fixture(params=[">", "<"])
    def server_client(self, universe, imdsinfo, port, request):
        server = InThreadIMDServer(universe.trajectory)
        imdsinfo.endianness = request.param
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", port, first_frame=False)
        client = IMDClient(
            f"localhost",
            port,
            universe.trajectory.n_atoms,
        )
        yield server, client
        client.stop()
        server.cleanup()

    def test_traj_unchanged(self, server_client, universe):
        server, client = server_client
        server.send_frames(0, 5)
        for i in range(5):
            imdf = client.get_imdframe()
            assert_allclose(universe.trajectory[i].time, imdf.time)
            assert_allclose(universe.trajectory[i].dt, imdf.dt)
            assert_allclose(universe.trajectory[i].data["step"], imdf.step)
            assert_allclose(universe.trajectory[i].positions, imdf.positions)
            assert_allclose(universe.trajectory[i].velocities, imdf.velocities)
            assert_allclose(universe.trajectory[i].forces, imdf.forces)
            assert_allclose(
                universe.trajectory[i].triclinic_dimensions, imdf.box
            )

    def test_pause_resume_continue(self, server_client_two_frame_buf):
        server, client = server_client_two_frame_buf
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

    def test_pause_resume_disconnect(self, server_client_two_frame_buf):
        """Client pauses because buffer is full, empties buffer and attempt to resume, but
        finds that simulation has already ended and raises EOF"""
        server, client = server_client_two_frame_buf
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

    def test_pause_resume_no_disconnect(self, server_client_two_frame_buf):
        """Client pauses because buffer is full, empties buffer and attempt to resume, but
        finds that simulation has already ended (but has not yet disconnected) and raises EOF
        """
        server, client = server_client_two_frame_buf
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
