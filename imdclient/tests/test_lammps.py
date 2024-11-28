import MDAnalysis as mda
import pytest
from pathlib import Path
import logging
from .base import IMDv3IntegrationTest
from .datafiles import LAMMPS_TOPOL, LAMMPS_IN_NST_1, LAMMPS_IN_NST_8

logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("lammps_test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


class TestIMDv3Lammps(IMDv3IntegrationTest):

    @pytest.fixture(params=[LAMMPS_IN_NST_1, LAMMPS_IN_NST_8])
    def inp(self, request):
        return request.param

    @pytest.fixture()
    def simulation_command(self, inp):
        return f"lmp < {Path(inp).name}"

    @pytest.fixture()
    def topol(self):
        return Path(LAMMPS_TOPOL).name

    @pytest.fixture()
    def traj(self):
        return "lammps_trj.h5md"

    @pytest.fixture()
    def input_files(self, inp):
        return [inp, LAMMPS_TOPOL]

    # @pytest.fixture()
    # def match_string(self):
    #     return "Waiting for IMD connection on port 8888"

    @pytest.fixture()
    def first_frame(self, inp):
        if inp == LAMMPS_IN_NST_1:
            return 1
        else:
            return 0

    # Not present in lammps-produced H5MD
    @pytest.fixture()
    def comp_dt(self):
        return False

    # This must wait until after imd stream has ended
    @pytest.fixture()
    def true_u(self, topol, traj, imd_u, tmp_path):
        u = mda.Universe(
            (tmp_path / topol),
            (tmp_path / traj),
            atom_style="id type x y z",
            convert_units=False,
        )
        yield u

    @pytest.fixture()
    def imd_u(self, docker_client, topol, tmp_path, port):
        u = mda.Universe(
            (tmp_path / topol),
            f"imd://localhost:{port}",
            atom_style="id type x y z",
        )
        with mda.Writer(
            (tmp_path / "imd.trr").as_posix(), u.trajectory.n_atoms
        ) as w:
            for ts in u.trajectory:
                w.write(u.atoms)
        yield mda.Universe(
            (tmp_path / topol),
            (tmp_path / "imd.trr"),
            atom_style="id type x y z",
        )
