import logging
from pathlib import Path
import re

import pytest
import MDAnalysis as mda

from .minimalreader import MinimalReader
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

    @pytest.fixture()
    def dt(self, inp):
        pattern = re.compile(r"^\s*timestep\s*(\S+)")
        with open(inp, "r") as file:
            for line in file:
                match = pattern.match(line)
                if match:
                    return float(match.group(1))
        raise ValueError(f"No dt found in {inp}")

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
        n_atoms = mda.Universe(
            tmp_path / topol,
            atom_style="id type x y z",
            convert_units=False,
        ).atoms.n_atoms
        u = MinimalReader(
            f"imd://localhost:{port}", n_atoms=n_atoms, process_stream=True
        )
        yield u
