import time
import logging
import select

logger = logging.getLogger("imdclient.IMDClient")


class timeit(object):
    """measure time spend in context

    :class:`timeit` is a context manager (to be used with the :keyword:`with`
    statement) that records the execution time for the enclosed context block
    in :attr:`elapsed`.

    Attributes
    ----------
    elapsed : float
        Time in seconds that elapsed between entering
        and exiting the context.

    Example
    -------
    Use as a context manager::

       with timeit() as total:
          # code to be timed

       print(total.elapsed, "seconds")

    See Also
    --------
    :func:`time.time`

    """

    def __enter__(self):
        self._start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.time()
        self.elapsed = end_time - self._start_time
        # always propagate exceptions forward
        return False


# NOTE: think of other edge cases as well- should be robust
def parse_host_port(filename):
    # Check if the format is correct
    parts = filename.split(":")
    if len(parts) == 2:
        host = parts[0]  # Hostname part
        try:
            port = int(parts[1])  # Convert the port part to an integer
            return (host, port)
        except ValueError:
            # Handle the case where the port is not a valid integer
            raise ValueError("Port must be an integer")
    else:
        # Handle the case where the format does not match "host:port"
        raise ValueError("Filename must be in the format 'host:port'")


def approximate_timestep_memsize(
    n_atoms, energies, dimensions, positions, velocities, forces
):
    total_size = 0

    if energies:
        total_size += 36

    if dimensions:
        # dimensions in the form (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        total_size += 24

    for dset in (positions, velocities, forces):
        if dset:
            total_size += n_atoms * 12

    return total_size


def read_into_buf(sock, buf):
    """Receives len(buf) bytes into buf from the socket sock"""
    view = memoryview(buf)
    total_received = 0
    while total_received < len(view):
        try:
            received = sock.recv_into(view[total_received:])
            if received == 0:
                logger.debug(
                    "read_into_buf excepting due to server closing connection"
                )
                raise ConnectionError
        except TimeoutError:
            logger.debug("read_into_buf excepting due to read timeout")
            raise TimeoutError
        except BlockingIOError:
            logger.debug(
                "read_into_buf excepting because socket timeout is 0 and no bytes are available to read"
            )
            raise BlockingIOError
        except OSError:
            logger.debug(
                "read_into_buf excepting because socket was closed elsewhere"
            )
            raise OSError
        except Exception as e:
            logger.debug(f"read_into_buf excepting due to: {e}")
            raise e
        total_received += received


def sock_contains_data(sock, timeout) -> bool:
    ready_to_read, ready_to_write, in_error = select.select(
        [sock], [], [], timeout
    )
    return sock in ready_to_read
