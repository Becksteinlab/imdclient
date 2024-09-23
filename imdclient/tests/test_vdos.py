import pytest
from imdclient.vdos import VDOS
import MDAnalysis as mda
from MDAnalysis.lib.log import ProgressBar


TOPOL = "/scratch/mheyden1/POPC/run-NVE.tpr"
TRAJ = "/scratch/mheyden1/POPC/run-NVE.trr"

# Source of truth

import numpy as np
from scipy.fft import fft


class vdos:
    def __init__(self, sel, nCorr):
        self.sel = sel
        self.nCorr = nCorr
        self.nRes = sel.residues.n_residues
        # initialize lists and arrays
        # Note: atMassLists are lists of numpy arrays
        #      => they are not a numpy arrays themselves
        self.atMassLists = (
            []
        )  # list of numpy arrays of atomic masses for each residue
        for res in self.sel.residues:
            # Note: np.newaxis is used to make the array 2D
            self.atMassLists.append(
                res.atoms.masses[:, np.newaxis].astype("float64")
            )
        # numpy array of residue masses from MDAnalysis
        self.resMassList = sel.residues.masses
        # numpy array buffers for COM position, COM velocity, angular momentum, and (sorted) moments of inertia
        self.COMposBuffer = np.zeros((nCorr, self.nRes, 3), dtype=np.float64)
        self.COMvelBuffer = np.zeros((nCorr, self.nRes, 3), dtype=np.float64)
        # time and frequency axes for correlation functions and VDoS
        self.tau = np.zeros(nCorr, dtype=np.float64)
        self.wavenumber = np.zeros(nCorr, dtype=np.float64)
        # numpy arrays for correlation functions and VDoS
        # -> rigid body translation
        self.trCorr = np.zeros((nCorr, self.nRes), dtype=np.float64)
        self.trVDoS = np.zeros((nCorr, self.nRes), dtype=np.float64)
        # initialize counter for normalization
        self.corrCnt = np.zeros(self.nRes, dtype=int)
        # flag for normalization
        self.normalized = 0

    def processStep(self, tStep, time):
        """process a single time step of the simulation"""
        idx = tStep % self.nCorr
        if tStep < self.nCorr:
            self.tau[tStep] = time
        pos = []
        vel = []
        r = 0
        # compute center of mass position and velocity
        for res in self.sel.residues:
            pos.append(res.atoms.positions.astype("float64"))
            vel.append(res.atoms.velocities.astype("float64"))
            # self.COMposBuffer[idx,r] = res.atoms.center_of_mass() # => too slow & no equivalent for velocity
            self.COMposBuffer[idx, r] = (
                np.sum(self.atMassLists[r] * pos[-1], axis=0)
                / self.resMassList[r]
            )
            self.COMvelBuffer[idx, r] = (
                np.sum(self.atMassLists[r] * vel[-1], axis=0)
                / self.resMassList[r]
            )
            r += 1
        # if sufficient data is available in buffers, compute correlation functions
        if tStep >= self.nCorr - 1:
            self.calcCorr(tStep + 1)

    def calcCorr(self, start):
        """compute correlation functions for all data in buffers"""
        # compute time correlation function for COM translation (for each residue)
        for i in range(self.nCorr):
            j = start % self.nCorr
            k = (j + i) % self.nCorr
            # for each residue
            self.trCorr[i] += np.sum(
                self.COMvelBuffer[j] * self.COMvelBuffer[k], axis=1
            )
        self.corrCnt += 1

    def normalize(self):
        """normalize correlation functions by number of data points"""
        if self.normalized == 0:
            for i in range(self.nRes):
                self.trCorr[:, i] *= self.resMassList[i] / self.corrCnt[i]
            self.normalized = 1

    def calcVDOS(self):
        """compute vibrational density of states from time correlation functions"""
        period = (self.tau[1] - self.tau[0]) * (2 * self.nCorr - 1)
        wn0 = (1.0 / period) * 33.35641
        self.wavenumber = np.arange(0, self.nCorr) * wn0
        tmp1 = np.zeros(2 * self.nCorr - 1, dtype=np.float64)
        for i in range(self.nRes):
            for j in range(self.nCorr):
                tmp1[j] = self.trCorr[j][i]
            for j in range(1, self.nCorr):
                k = 2 * self.nCorr - j - 1
                tmp1[k] = tmp1[j]
            tmp1 = fft(tmp1)
            for j in range(self.nCorr):
                self.trVDoS[j][i] = tmp1[j].real


@pytest.fixture
def true_vdos():
    u = mda.Universe(TOPOL, TRAJ)
    sel = u.select_atoms("resname POPC")
    true_vdos = vdos(sel, 200)
    tStep = 0
    for ts in ProgressBar(u.trajectory[:400]):
        true_vdos.processStep(tStep, ts.time)
        tStep += 1
    true_vdos.normalize()
    true_vdos.calcVDOS()
    yield (
        true_vdos.trCorr,
        true_vdos.trVDoS,
        true_vdos.wavenumber,
        true_vdos.tau,
    )


@pytest.fixture
def mda_vdos():
    u = mda.Universe(TOPOL, TRAJ)
    sel = u.select_atoms("resname POPC")
    vdos = VDOS(u.trajectory, sel, nCorr=200).run(stop=400)
    yield (
        vdos.results.trCorr,
        vdos.results.trVDoS,
        vdos.results.wavenumber,
        vdos.results.tau,
    )


def test_compare_vdos(true_vdos, mda_vdos):
    for true, mda in zip(true_vdos, mda_vdos):
        np.testing.assert_allclose(true, mda, rtol=1e-3)
