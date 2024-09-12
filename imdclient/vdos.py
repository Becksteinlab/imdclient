import MDAnalysis as mda
from scipy.fft import rfft, fft
import numpy as np
from MDAnalysis.analysis.base import (
    AnalysisBase,
)
from MDAnalysis.lib.log import ProgressBar

import logging

logger = logging.getLogger(__name__)


class StreamFriendlyAnalysisBase(AnalysisBase):
    _analysis_algorithm_is_parallelizable = False

    def __init__(self, trajectory, verbose=False, **kwargs):
        super().__init__(trajectory, verbose=verbose, **kwargs)
        if hasattr(trajectory, "one_pass") and trajectory.one_pass:
            self._one_pass = True
        else:
            self._one_pass = False

    def run(
        self,
        start=None,
        stop=None,
        step=None,
        frames=None,
        verbose=None,
        *,
        progressbar_kwargs={},
    ):
        """Perform the calculation

        Parameters
        ----------
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        frames : array_like, optional
            array of integers or booleans to slice trajectory; `frames` can
            only be used *instead* of `start`, `stop`, and `step`. Setting
            *both* `frames` and at least one of `start`, `stop`, `step` to a
            non-default value will raise a :exc:`ValueError`.

            .. versionadded:: 2.2.0

        verbose : bool, optional
            Turn on verbosity

        progressbar_kwargs : dict, optional
            ProgressBar keywords with custom parameters regarding progress bar position, etc;
            see :class:`MDAnalysis.lib.log.ProgressBar` for full list.


        .. versionchanged:: 2.2.0
            Added ability to analyze arbitrary frames by passing a list of
            frame indices in the `frames` keyword argument.

        .. versionchanged:: 2.5.0
            Add `progressbar_kwargs` parameter,
            allowing to modify description, position etc of tqdm progressbars
        """
        logger.info("Choosing frames to analyze")
        # if verbose unchanged, use class default
        verbose = (
            getattr(self, "_verbose", False) if verbose is None else verbose
        )

        if self._one_pass:
            self._sliced_trajectory = self._trajectory
            self.frames = np.zeros(1000, dtype=int)
            self.times = np.zeros(1000)
        else:
            self._setup_frames(
                self._trajectory,
                start=start,
                stop=stop,
                step=step,
                frames=frames,
            )

        logger.info("Starting preparation")
        self._prepare()

        if self._one_pass:
            logger.info(
                "Starting analysis loop over trajectory frames",
            )
        else:
            logger.info(
                "Starting analysis loop over %d trajectory frames",
                self.n_frames,
            )

        for i, ts in enumerate(
            ProgressBar(
                self._sliced_trajectory, verbose=verbose, **progressbar_kwargs
            )
        ):
            self._frame_index = i
            self._ts = ts
            if self._one_pass:
                self._alloc_frame_and_time_arrays()
            self.frames[i] = ts.frame
            self.times[i] = ts.time
            self._single_frame()
        logger.info("Finishing up")

        if self._one_pass:
            self.frames = self.frames[: self._frame_index + 1]
            self.times = self.times[: self._frame_index + 1]
            self.n_frames = self._frame_index + 1

        self._conclude()

        return self

    def _alloc_frame_and_time_arrays(self):
        if self._frame_index == len(self.frames) - 1:
            self.frames = np.resize(self.frames, 1000 + len(self.frames))
            self.times = np.resize(self.times, 1000 + len(self.times))


class VDOS(StreamFriendlyAnalysisBase):

    def __init__(
        self, trajectory, selection, nCorr=100, verbose=False, **kwargs
    ):
        super().__init__(trajectory, verbose=False, **kwargs)
        self.sel = selection
        self.nCorr = nCorr
        self.nRes = selection.residues.n_residues

    def _prepare(self):
        # numpy array buffers for COM position, COM velocity, angular momentum, and (sorted) moments of inertia
        self.COMposBuffer = np.zeros(
            (self.nCorr, self.nRes, 3), dtype=np.float64
        )
        self.COMvelBuffer = np.zeros(
            (self.nCorr, self.nRes, 3), dtype=np.float64
        )

        # time and frequency axes for correlation functions and VDoS
        self.results["tau"] = np.zeros(self.nCorr, dtype=np.float64)
        self.results["wavenumber"] = np.zeros(self.nCorr, dtype=np.float64)
        # numpy arrays for correlation functions and VDoS
        # -> rigid body translation
        self.results["trCorr"] = np.zeros(
            (self.nCorr, self.nRes), dtype=np.float64
        )
        self.results["trVDoS"] = np.zeros(
            (self.nCorr, self.nRes), dtype=np.float64
        )
        # initialize counter for normalization
        self.corrCnt = 0

        # # If residues can have different num atoms, use list
        # self.atMassLists = []
        # for res in self.sel.residues:
        #     self.atMassLists.append(
        #         res.atoms.masses[:, np.newaxis].astype("float64")
        #     )

        self.atMassLists = np.empty(
            (len(self.sel.residues), len(self.sel.residues[0].atoms), 1),
            dtype=np.float64,
        )
        for i, res in enumerate(self.sel.residues):
            self.atMassLists[i] = res.atoms.masses[:, np.newaxis].astype(
                "float64"
            )

    def _single_frame(self):
        idx = self._ts.frame % self.nCorr

        if self._ts.frame < self.nCorr:
            self.results.tau[idx] = self._ts.time

        atMassLists = self.atMassLists
        sel = self.sel
        masses = sel.residues.masses
        COMvelBuffer = self.COMvelBuffer

        # # compute center of mass position and velocity
        for i in range(self.nRes):
            vel = sel.residues[i].atoms.velocities.astype("float64")
            COMvelBuffer[idx, i] = (
                np.sum(atMassLists[i] * vel, axis=0) / masses[i]
            )

        # if sufficient data is available in buffers, compute correlation functions
        if self._ts.frame >= self.nCorr - 1:
            self._calcCorr(self._ts.frame + 1)

    def _calcCorr(self, start):
        """compute correlation functions for all data in buffers"""
        # compute time correlation function for COM translation (for each residue)
        # can be parallelized
        trCorr = self.results.trCorr
        COMvelBuffer = self.COMvelBuffer
        nCorr = self.nCorr

        for i in range(nCorr):
            j = start % nCorr
            k = (j + i) % nCorr
            trCorr[i] += np.sum(COMvelBuffer[j] * COMvelBuffer[k], axis=1)
        self.corrCnt += 1

    def _conclude(self):
        """normalize correlation functions by number of data points
        and calculate compute vibrational density of states from time correlation functions
        """
        # Normalization
        # multiply by mass to convert velocity -> momentum
        self.results.trCorr[:] *= self.sel.residues.masses[:] / self.corrCnt
        ##  Calculate VDoS
        period = (self.results.tau[1] - self.results.tau[0]) * (
            2 * self.nCorr - 1
        )
        wn0 = (1.0 / period) * 33.35641
        self.results.wavenumber = np.arange(0, self.nCorr) * wn0
        tmp1 = np.zeros((2 * self.nCorr - 1, self.nRes), dtype=np.float64)
        tmp1[: self.nCorr] = self.results.trCorr[:]
        tmp1[self.nCorr :] = tmp1[1 : self.nCorr][::-1]
        self.results.trVDoS = rfft(tmp1, axis=0)
