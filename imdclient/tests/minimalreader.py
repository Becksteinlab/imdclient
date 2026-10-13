import logging
import copy

import numpy as np
from MDAnalysis.coordinates import core
from MDAnalysis.lib import distances

from imdclient.IMDClient import IMDClient
from imdclient.utils import parse_host_port

logger = logging.getLogger("imdclient.IMDClient")


class MinimalReader:
    """
    Minimal reader for testing purposes

    Parameters
    ----------
    filename : str
        a string of the form "host:port" where host is the hostname
        or IP address of the listening MD engine server and port
        is the port number.
    n_atoms : int
        number of atoms in the system. defaults to number of atoms
        in the topology. don't set this unless you know what you're doing
    process_stream : bool (optional)
        if True, the reader will process the stream of frames and
        store them in the `trajectory` attribute. defaults to False.
    kwargs : dict (optional)
        keyword arguments passed to the constructed :class:`IMDClient`
    """

    def __init__(self, filename, n_atoms, process_stream=False, **kwargs):

        self.imd_frame = None

        # a trajectory of imd frames
        self.trajectory = []

        self.n_atoms = n_atoms

        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = IMDClient(host, port, n_atoms, **kwargs)

        self._frame = -1

        if process_stream:
            self._process_stream()

    @property
    def imdsinfo(self):
        return self._imdclient.get_imdsessioninfo()

    def _read_next_frame(self):
        try:
            imd_frame = self._imdclient.get_imdframe()
        except EOFError:
            raise

        self._frame += 1
        self.imd_frame = imd_frame

        # Modify the box dimensions to be triclinic
        self._modify_box_dimensions()

        logger.debug(f"MinimalReader: Loaded frame {self._frame}")

        return self.imd_frame

    def _modify_box_dimensions(self):
        if self.imd_frame.box is None:
            self.imd_frame.dimensions = None
            return
        self.imd_frame.dimensions = core.triclinic_box(*self.imd_frame.box)

    def _process_stream(self):
        # Process the stream of frames
        while True:
            try:
                self.trajectory.append(copy.deepcopy(self._read_next_frame()))
                # `.copy()` might not be required but adding it to cover any edge cases where a refernce gets passed
                logger.debug(
                    f"MinimalReader: Added frame {self._frame} to trajectory"
                )
            except EOFError:
                break

    def _wrap_frame(self, frame, box):
        if box is None:
            return
        frame.positions = distances.apply_PBC(
            np.asarray(frame.positions, dtype=np.float32),
            box,
        )

    def _wrap_trajectory(self, true_u=None, first_frame=0):
        """Wrap each stored IMD frame into the primary box (per atom)."""
        if true_u is None:
            for frame in self.trajectory:
                self._wrap_frame(frame, frame.dimensions)
            return

        for i in range(first_frame, len(true_u.trajectory)):
            self._wrap_frame(
                self.trajectory[i - first_frame],
                true_u.trajectory[i].dimensions,
            )

    def close(self):
        """Gracefully shut down the reader. Stops the producer thread."""
        logger.debug("MinimalReader: close() called")
        if self._imdclient is not None:
            self._imdclient.stop()
