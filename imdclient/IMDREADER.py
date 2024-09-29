"""

Example: Streaming an IMD v2 trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To stream a trajectory from GROMACS or another simulation engine that supports 
IMD v2, ensure that the simulation engine is running and waiting for an IMD connection.

For example, in GROMACS, you can use ``gmx mdrun`` with the ``-imdwait`` flag
to ensure that GROMACS will wait for a client before starting the simulation.
In GROMACS, you will know that the simulation is ready and waiting for the
MDAnalysis IMDReader client when this line is printed to the terminal:

.. code-block:: none

    IMD: Will wait until I have a connection and IMD_GO orders.

Once the simulation is ready for a client connection, setup your :class:`Universe`
like this: ::

    import MDAnalysis as mda
    # Pass host and port of the listening GROMACACS simulation
    # server as the trajectory argument
    u = mda.Universe("topology.tpr", "localhost:8888")

For more information on IMD v2 as it is implemented in GROMACS, see the imd command line
arguments described `here <https://manual.gromacs.org/documentation/5.1/onlinehelp/gmx-mdrun.html>`_,
the :ref:`v2_spec` writeup, and this `example simulation  <https://www.mpinat.mpg.de/grubmueller/interactivemd>`_. Note that 
setting the following line in your ``.mdp`` file allows you to change which subset of the
simulation is sent to the client via IMD v2:

.. code-block:: none

    ; Group to display and/or manipulate in interactive MD session
    IMD-group                = System

Classes
^^^^^^^

.. autoclass:: IMDReader
   :members:
   :inherited-members:

"""

from MDAnalysis.coordinates import core
from MDAnalysis.lib.util import store_init_arguments

# NOTE: changeme
from .IMDClient import IMDClient
from .utils import *
import numpy as np
import logging

from typing import Optional

from .streambase import StreamReaderBase

logger = logging.getLogger("imdclient.IMDClient")


class IMDReader(StreamReaderBase):
    """
    Reader for IMD protocol packets.
    """

    format = "IMD"
    one_pass = True

    @store_init_arguments
    def __init__(
        self,
        filename,
        convert_units=True,
        n_atoms=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        filename : a string of the form "host:port" where host is the hostname
            or IP address of the listening GROMACS server and port
            is the port number.
        n_atoms : int (optional)
            number of atoms in the system. defaults to number of atoms
            in the topology. don't set this unless you know what you're doing.
        """

        super(IMDReader, self).__init__(filename, **kwargs)

        logger.debug("IMDReader initializing")

        if n_atoms is None:
            raise ValueError(
                "IMDReader: n_atoms must be specified"
            )
        self.n_atoms = n_atoms

        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = IMDClient(
            host, port, n_atoms, **kwargs
        )

        imdsinfo = self._imdclient.get_imdsessioninfo()
        # NOTE: after testing phase, fail out on IMDv2

        self.ts = self._Timestep(
            self.n_atoms,
            positions=imdsinfo.positions,
            velocities=imdsinfo.velocities,
            forces=imdsinfo.forces,
            **self._ts_kwargs,
        )

        self._frame = -1

        try:
            self._read_next_timestep()
        except StopIteration:
            raise RuntimeError(
                "IMDReader: No data found in stream"
            )

    def _read_frame(self, frame):

        try:
            imdf = self._imdclient.get_imdframe()
        except EOFError:
            # Not strictly necessary, but for clarity
            raise StopIteration

        self._frame = frame
        self._load_imdframe_into_ts(imdf)

        logger.debug(
            f"IMDReader: Loaded frame {self._frame}"
        )
        return self.ts

    def _load_imdframe_into_ts(self, imdf):
        self.ts.frame = self._frame
        if imdf.time is not None:
            self.ts.time = imdf.time
            # NOTE: timestep.pyx "dt" method is suspicious bc it uses "new" keyword for a float
            self.ts.data["dt"] = imdf.dt
            self.ts.data["step"] = imdf.step
        if imdf.energies is not None:
            self.ts.data.update(imdf.energies)
        if imdf.box is not None:
            self.ts.dimensions = core.triclinic_box(
                *imdf.box
            )
        if imdf.positions is not None:
            # must call copy because reference is expected to reset
            # see 'test_frame_collect_all_same' in MDAnalysisTests.coordinates.base
            self.ts.positions = imdf.positions
        if imdf.velocities is not None:
            self.ts.velocities = imdf.velocities
        if imdf.forces is not None:
            self.ts.forces = imdf.forces

    @staticmethod
    def _format_hint(thing):
        try:
            parse_host_port(thing)
        except:
            return False
        return True

    def close(self):
        """Gracefully shut down the reader. Stops the producer thread."""
        logger.debug("IMDReader close() called")
        self._imdclient.stop()
        # NOTE: removeme after testing
        logger.debug("IMDReader shut down gracefully.")
