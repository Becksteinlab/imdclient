import logging
import select
import time
from .IMDClient import IMDClient

logger = logging.getLogger("imdclient.minimalReader")

class minimalReader:
    """
    Minimal reader for testing purposes
    
    Parameters
    ----------
    filename : str
        a string of the form "host:port" where host is the hostname
        or IP address of the listening GROMACS server and port
        is the port number.
        n_atoms : int (optional)
        number of atoms in the system. defaults to number of atoms
        in the topology. don't set this unless you know what you're doing.
        kwargs : dict (optional)
        keyword arguments passed to the constructed :class:`IMDClient`
        """
    def __init__(self, filename, n_atoms=None, **kwargs):
        self._init_scope = True
        
        self._first_df = None

        self.imdf = None

        self.n_atoms = n_atoms
        
        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = IMDClient(host, port, n_atoms, **kwargs)

        imdsinfo = self._imdclient.get_imdsessioninfo()

        self._frame = -1

    def _read_frame_into_df(self, frame):

        try:
            imdf = self._imdclient.get_imdframe()
        except EOFError as e:
            raise e

        self._frame = frame
        self.imdf = imdf

        logger.debug(f"minimalReader: Loaded frame {self._frame}")
        return imdf

    def _read_next_timestep(self):
        # No rewinding- to both load the first frame after  __init__
        # and access it again during iteration, we need to store first ts in mem
        if not self._init_scope and self._frame == -1:
            self._frame += 1
            # can't simply return the same ts again- transformations would be applied twice
            # instead, return the pre-transformed copy
            return self._first_df

        self._read_frame_into_df(self._frame + 1)

        if self._init_scope:
            self._first_df = self.imdf.copy()
            self._init_scope = False

        return imdf

def parse_host_port(filename):
    if not filename.startswith("imd://"):
        raise ValueError("minimalReader: URL must be in the format 'imd://host:port'")
    
    # Check if the format is correct
    parts = filename.split("imd://")[1].split(":")
    if len(parts) == 2:
        host = parts[0] 
        try:
            port = int(parts[1])
            return (host, port)
        except ValueError:
            raise ValueError("minimalReader: Port must be an integer")
    else:
        raise ValueError("minimalReader: URL must be in the format 'imd://host:port'")
