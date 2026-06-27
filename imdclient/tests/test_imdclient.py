"""Test for IMDClient functionality"""

import logging

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (
    COORDINATES_H5MD,
    COORDINATES_TOPOLOGY,
)
from numpy.testing import assert_allclose
import pytest

from imdclient.IMDClient import IMDClient, imdframe_memsize
from imdclient.IMDProtocol import (
    IMDHeaderType,
    create_energy_bytes,
    create_header_bytes,
)
from .server import InThreadIMDServer
from .utils import (
    create_default_imdsinfo_v2,
    create_default_imdsinfo_v3,
)


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


class IMDClientTest:
    @pytest.fixture
    def universe(self):
        return mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)

    @pytest.fixture
    def imdsinfo(self):
        return create_default_imdsinfo_v2()

    @pytest.fixture
    def server_client_two_frame_buf(self, universe, imdsinfo):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=False)
        client = IMDClient(
            "localhost",
            server.port,
            universe.trajectory.n_atoms,
            buffer_size=imdframe_memsize(universe.trajectory.n_atoms, imdsinfo)
            * 2,
        )
        server.join_accept_thread()
        yield server, client
        client.stop()
        server.cleanup()

    @pytest.fixture(params=[">", "<"])
    def server_client(self, universe, imdsinfo, request):
        server = InThreadIMDServer(universe.trajectory)
        imdsinfo.endianness = request.param
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=False)
        client = IMDClient("localhost", server.port, universe.atoms.n_atoms)
        server.join_accept_thread()
        yield server, client
        client.stop()
        server.cleanup()

    @pytest.fixture
    def server_client_incorrect_atoms(self, universe, imdsinfo):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=False)
        client = IMDClient("localhost", server.port, universe.atoms.n_atoms + 1)
        server.join_accept_thread()
        yield server, client
        client.stop()
        server.cleanup()

    def test_traj_unchanged(self, server_client, imdsinfo, universe):
        server, client = server_client
        server.send_frames(0, 5)
        for i in range(5):
            imdf = client.get_imdframe()
            if imdsinfo.time:
                assert_allclose(universe.trajectory[i].time, imdf.time)
                assert_allclose(universe.trajectory[i].dt, imdf.dt)
                assert_allclose(universe.trajectory[i].data["step"], imdf.step)
            if imdsinfo.box:
                assert_allclose(
                    universe.trajectory[i].triclinic_dimensions, imdf.box
                )
            if imdsinfo.positions:
                assert_allclose(
                    universe.trajectory[i].positions, imdf.positions
                )
            if imdsinfo.velocities:
                assert_allclose(
                    universe.trajectory[i].velocities, imdf.velocities
                )
            if imdsinfo.forces:
                assert_allclose(universe.trajectory[i].forces, imdf.forces)

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


class TestIMDClientV3(IMDClientTest):
    @pytest.fixture
    def imdsinfo(self):
        return create_default_imdsinfo_v3()

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
    def test_continue_after_disconnect(self, universe, imdsinfo, cont):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=False)
        client = IMDClient(
            "localhost",
            server.port,
            universe.trajectory.n_atoms,
            continue_after_disconnect=cont,
        )
        server.join_accept_thread()
        server.expect_packet(
            IMDHeaderType.IMD_WAIT, expected_length=int(not cont)
        )

    def test_trate_not_sent_for_v3(self, universe, imdsinfo):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=False)
        client = IMDClient(
            "localhost",
            server.port,
            universe.atoms.n_atoms,
            transmission_rate=8,
        )
        server.join_accept_thread()
        server.expect_no_packet()
        client.stop()
        server.cleanup()


class TestIMDClientV2(IMDClientTest):
    @pytest.mark.parametrize("rate", [1, 8])
    def test_trate_sent_after_go_for_v2(self, universe, imdsinfo, rate):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", first_frame=False)
        client = IMDClient(
            "localhost",
            server.port,
            universe.atoms.n_atoms,
            transmission_rate=rate,
        )
        server.join_accept_thread()
        server.expect_packet(IMDHeaderType.IMD_TRATE, expected_length=rate)
        client.stop()
        server.cleanup()

    def test_pause_pause_continue(self, server_client_two_frame_buf):
        server, client = server_client_two_frame_buf
        server.send_frames(0, 2)
        # Client's buffer is filled; client should send pause.
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        # Empty buffer.
        client.get_imdframe()
        # Only the second call actually frees buffer memory.
        client.get_imdframe()
        # IMDv2 uses IMD_PAUSE for both pause and unpause.
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        server.send_frame(1)
        client.get_imdframe()

    def test_pause_pause_disconnect(self, server_client_two_frame_buf):
        """Client pauses because buffer is full, empties buffer and attempts to
        unpause with a second IMD_PAUSE, but simulation has already ended and
        raises EOF."""
        server, client = server_client_two_frame_buf
        server.send_frames(0, 2)
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        client.get_imdframe()
        client.get_imdframe()
        # IMDv2 uses IMD_PAUSE for both pause and unpause.
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        # Simulation is over; client should raise EOF.
        server.disconnect()
        with pytest.raises(EOFError):
            client.get_imdframe()

    def test_pause_pause_no_disconnect(self, server_client_two_frame_buf):
        """Client pauses because buffer is full, empties buffer and attempts to
        unpause with a second IMD_PAUSE, but simulation has already ended (without
        disconnecting) and raises EOF."""
        server, client = server_client_two_frame_buf
        server.send_frames(0, 2)
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        client.get_imdframe()
        client.get_imdframe()
        # IMDv2 uses IMD_PAUSE for both pause and unpause.
        server.expect_packet(IMDHeaderType.IMD_PAUSE)
        # Simulation is over; client should raise EOF.
        with pytest.raises(EOFError):
            client.get_imdframe()
        # Server should receive disconnect from client (though it doesn't have to do anything).
        server.expect_packet(IMDHeaderType.IMD_DISCONNECT)

    def test_reads_multiple_leading_energies_before_coords(
        self, server_client, caplog
    ):
        server, client = server_client

        with caplog.at_level(logging.WARNING, logger="imdclient.IMDClient"):
            endianness = client.get_imdsessioninfo().endianness
            # Send >3 leading energy packets before coordinates.
            for i in range(4):
                energy_header = create_header_bytes(
                    IMDHeaderType.IMD_ENERGIES, 1
                )
                energies = create_energy_bytes(
                    i,
                    i + 1,
                    i + 2,
                    i + 3,
                    i + 4,
                    i + 5,
                    i + 6,
                    i + 7,
                    i + 8,
                    i + 9,
                    endianness,
                )
                server.conn.sendall(energy_header + energies)

            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, server.traj.n_atoms
            )
            pos = universe_frame0 = server.traj[0].positions
            pos_bytes = pos.astype(f"{endianness}f", copy=False).tobytes()
            server.conn.sendall(pos_header + pos_bytes)

            imdf = client.get_imdframe()

        assert imdf.energies is not None
        assert imdf.energies["step"] == 3
        assert_allclose(universe_frame0, imdf.positions)
        assert any(
            "Received 4 leading IMDv2 energy packets" in rec.message
            for rec in caplog.records
        )

    def test_reads_no_leading_energy_packets(self, server_client, caplog):
        server, client = server_client

        with caplog.at_level(logging.WARNING, logger="imdclient.IMDClient"):
            endianness = client.get_imdsessioninfo().endianness
            # Send coords only, no leading energy packet.
            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, server.traj.n_atoms
            )
            pos = server.traj[0].positions
            pos_bytes = pos.astype(f"{endianness}f", copy=False).tobytes()
            server.conn.sendall(pos_header + pos_bytes)

            imdf = client.get_imdframe()

        assert imdf.energies is None
        assert_allclose(pos, imdf.positions)
        assert any(
            "Received 0 leading IMDv2 energy packets" in rec.message
            for rec in caplog.records
        )

    def test_uses_previous_energies_when_frame_has_coords_only(
        self, server_client
    ):
        server, client = server_client

        # First frame with energies sets the cache.
        server.send_frame(0)
        first = client.get_imdframe()
        assert first.energies is not None
        prev_step = first.energies["step"]

        # Second frame: send only coords; should reuse previous energies.
        endianness = client.get_imdsessioninfo().endianness
        pos_header = create_header_bytes(
            IMDHeaderType.IMD_FCOORDS, server.traj.n_atoms
        )
        pos = server.traj[1].positions
        pos_bytes = pos.astype(f"{endianness}f", copy=False).tobytes()
        server.conn.sendall(pos_header + pos_bytes)

        second = client.get_imdframe()
        assert second.energies is not None
        assert second.energies["step"] == prev_step
        assert_allclose(pos, second.positions)

    @pytest.mark.parametrize(
        "packet_type",
        [
            IMDHeaderType.IMD_TIME,
            IMDHeaderType.IMD_BOX,
            IMDHeaderType.IMD_VELOCITIES,
            IMDHeaderType.IMD_FORCES,
        ],
    )
    def test_unexpected_packet_type_in_v2(self, server_client, packet_type):
        """IMDv2 rejects v3-only packet types received before coordinates."""
        server, client = server_client

        server.conn.sendall(create_header_bytes(packet_type, 1))

        with pytest.raises(EOFError) as exc_info:
            client.get_imdframe()

        assert f"Unexpected packet type {packet_type.name}" in str(
            exc_info.value
        )


class TestIMDClientV3ContextManager(IMDClientTest):
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
            "localhost", server.port, universe.trajectory.n_atoms
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
