import logging
import select
import time
from imdclient.IMDClient import IMDClient
from .utils import parse_host_port

logger = logging.getLogger("imdclient.IMDClient")


class minimalReader:
    """
    Minimal reader for testing purposes

    Parameters
    ----------
    filename : str
        a string of the form "host:port" where host is the hostname
        or IP address of the listening MD engine server and port
        is the port number.
    n_atoms : int (optional)
        number of atoms in the system. defaults to number of atoms
        in the topology. don't set this unless you know what you're doing.
    kwargs : dict (optional)
        keyword arguments passed to the constructed :class:`IMDClient`
    """

    def __init__(self, filename, **kwargs):

        self.imd_frame = None

        # a trajectory of imd frames
        self.trajectory = []

        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = IMDClient(host, port, n_atoms, **kwargs)

        imdsinfo = self._imdclient.get_imdsessioninfo()

        self._frame = -1

        self._process_stream()

    def _read_next_frame(self):
        try:
            imd_frame = self._imdclient.get_imdframe()
        except EOFError:
            raise

        self._frame += 1
        self.imd_frame = imd_frame

        logger.debug(f"minimalReader: Loaded frame {self._frame}")

        return imd_frame

    def _process_stream(self):
        # Process the stream of frames
        while True:
            try:
                self.trajectory.append(self._read_next_frame().copy())
                logger.debug(
                    f"minimalReader: Added frame {self._frame} to trajectory"
                )
            except EOFError:
                break
