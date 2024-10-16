import imdclient
from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_H5MD,
)
import MDAnalysis as mda
from .utils import (
    get_free_port,
    create_default_imdsinfo_v3,
)
from .server import InThreadIMDServer
from MDAnalysis.analysis.rms import RMSF
from numpy.testing import (
    assert_almost_equal,
)
import numpy as np
import pytest
from imdclient.IMDREADER import IMDReader


class TestStreamAnalysis:

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
    def imd_universe(self, universe, imdsinfo, port):
        server = InThreadIMDServer(universe.trajectory)
        server.set_imdsessioninfo(imdsinfo)
        server.handshake_sequence("localhost", port, first_frame=True)

        imd_universe = mda.Universe(COORDINATES_TOPOLOGY, f"localhost:{port}")
        server.send_frames(1, 5)

        yield imd_universe
        server.cleanup()

    def test_rmsf(self, imd_universe, universe):
        imd_rmsf = RMSF(imd_universe.atoms).run()
        rmsf = RMSF(universe.atoms).run()

        assert_almost_equal(imd_rmsf.results.rmsf, rmsf.results.rmsf)

    def test_stack_rmsf(self, imd_universe, universe):
        r1 = RMSF(imd_universe.atoms)
        r2 = RMSF(imd_universe.atoms)
        imdclient.StackableAnalysis(imd_universe.trajectory, [r1, r2]).run()

        rmsf = RMSF(universe.atoms).run()

        assert_almost_equal(r1.results.rmsf, rmsf.results.rmsf)
        assert_almost_equal(r2.results.rmsf, rmsf.results.rmsf)
