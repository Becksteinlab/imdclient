"""Test for IMDClient functionality"""

import logging
import time

import pytest
from numpy.testing import (
    assert_allclose,
)
import MDAnalysis as mda
from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_H5MD,
)

from imdclient.IMDClient import imdframe_memsize, IMDClient
from imdclient.IMDProtocol import IMDHeaderType
from .utils import (
    create_default_imdsinfo_v3,
)
from .server import InThreadIMDServer


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
    def universe(self):
        return mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)

    @pytest.fixture
    def imdsinfo(self):
        return create_default_imdsinfo_v3()

    @pytest.fixture
    def server_client(self, universe, imdsinfo):
        created = []

        def _server_client(endianness=None, **client_kwargs):
            server = InThreadIMDServer(universe.trajectory)
            if endianness is not None:
                imdsinfo.endianness = endianness
            server.set_imdsessioninfo(imdsinfo)

            n_atoms = client_kwargs.pop("n_atoms", universe.atoms.n_atoms)
            server.handshake_sequence("localhost", first_frame=False)
            client = IMDClient(
                "localhost",
                server.port,
                n_atoms,
                **client_kwargs,
            )
            server.join_accept_thread()
            created.append((server, client))
            return server, client

        yield _server_client

        # Teardown: stop clients and cleanup servers
        for server, client in created:
            try:
                client.stop()
            except Exception:
                pass
            try:
                server.cleanup()
            except Exception:
                pass

    @pytest.fixture
    def server_client_two_frame_buf(self, server_client, universe, imdsinfo):
        # Calculate the buffer size
        buffer_size = (
            imdframe_memsize(universe.trajectory.n_atoms, imdsinfo) * 2
        )
        timeout = 5 # to speed up no disconnect test
        server, client = server_client(buffer_size=buffer_size, timeout=timeout)
        yield server, client

    @pytest.fixture(params=["<", ">"])
    def server_client_endianness(self, server_client, request):
        server, client = server_client(endianness=request.param)
        yield server, client

    @pytest.fixture
    def server_client_incorrect_atoms(self, server_client, universe):
        server, client = server_client(n_atoms=universe.trajectory.n_atoms + 1)
        yield server, client

    def test_traj_unchanged(self, server_client_endianness, universe):
        server, client = server_client_endianness
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

    @pytest.mark.parametrize("cont", [True, False])
    def test_continue_after_disconnect(self, server_client, cont):
        server, client = server_client(continue_after_disconnect=cont)
        server.expect_packet(
            IMDHeaderType.IMD_WAIT, expected_length=(int)(not cont)
        )

    def test_timeout_warning_low_value(self, server_client, caplog):
        """Test that warning is issued for timeout values <= 1 second"""
        with caplog.at_level(logging.WARNING):
            server, client = server_client(timeout=1)

        # Check that warning was logged
        assert any(
            "timeout value of 1 second(s) is very low" in record.message
            for record in caplog.records
        )

    @pytest.mark.parametrize("timeout_val", [2, 10])
    def test_timeout_within_limit(self, server_client, universe, timeout_val):
        """Test that timeout does not trigger when server responds within timeout period"""
        server, client = server_client(timeout=timeout_val)

        # Sleep for less than timeout before sending frames
        time.sleep(timeout_val - 1)
        server.send_frame(0)

        # Should successfully receive frame without timeout
        imdf = client.get_imdframe()
        assert_allclose(universe.trajectory[0].positions, imdf.positions)

    @pytest.mark.parametrize("timeout_val", [2, 10])
    def test_timeout_when_exceeded(self, server_client, timeout_val):
        """Test that timeout triggers EOFError when server doesn't respond within timeout period"""
        server, client = server_client(timeout=timeout_val)

        # Sleep for longer than timeout without sending any frames
        time.sleep(timeout_val + 1)

        # Client should timeout and raise EOFError when trying to get first frame
        with pytest.raises(EOFError) as exc_info:
            client.get_imdframe()

        # Verify TimeoutError is somewhere in the exception chain
        exception_chain = []
        current = exc_info.value
        while current is not None:
            exception_chain.append(type(current))
            current = current.__cause__

        assert TimeoutError in exception_chain

    def test_incorrect_atom_count(
        self, server_client_incorrect_atoms, universe
    ):
        server, client = server_client_incorrect_atoms

        server.send_frame(0)

        with pytest.raises(EOFError) as exc_info:
            client.get_imdframe()

        error_msg = str(exc_info.value)
        assert (
            f"Expected n_atoms value {universe.atoms.n_atoms + 1}" in error_msg
        )
        assert f"got {universe.atoms.n_atoms}" in error_msg
        assert "Ensure you are using the correct topology file" in error_msg


class TestIMDClientV3ContextManager:
    @pytest.fixture
    def universe(self):
        return mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)

    @pytest.fixture
    def imdsinfo(self):
        return create_default_imdsinfo_v3()

    @pytest.fixture
    def server(self, universe, imdsinfo):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        yield server
        server.cleanup()

    def test_context_manager_traj_unchanged(self, server, universe):
        server.handshake_sequence("localhost", first_frame=False)

        i = 0
        with IMDClient(
            "localhost",
            server.port,
            universe.trajectory.n_atoms,
        ) as client:
            server.send_frames(0, 5)
            while i < 5:

                imdf = client.get_imdframe()
                assert_allclose(universe.trajectory[i].time, imdf.time)
                assert_allclose(universe.trajectory[i].dt, imdf.dt)
                assert_allclose(universe.trajectory[i].data["step"], imdf.step)
                assert_allclose(
                    universe.trajectory[i].positions, imdf.positions
                )
                assert_allclose(
                    universe.trajectory[i].velocities, imdf.velocities
                )
                assert_allclose(universe.trajectory[i].forces, imdf.forces)
                assert_allclose(
                    universe.trajectory[i].triclinic_dimensions, imdf.box
                )
                i += 1
        server.expect_packet(IMDHeaderType.IMD_DISCONNECT)
        assert i == 5
