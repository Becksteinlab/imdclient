from MDAnalysis.coordinates.base import (
    ReaderBase,
    FrameIteratorBase,
    FrameIteratorAll,
)
import numbers
import warnings


class StreamReaderBase(ReaderBase):

    def __init__(self, filename, convert_units=True, **kwargs):
        super(StreamReaderBase, self).__init__(
            filename, convert_units=convert_units, **kwargs
        )
        self._init_scope = True
        self._reopen_called = False
        self._first_ts = None

    def _read_next_timestep(self):
        # No rewinding- to both load the first frame after  __init__
        # and access it again during iteration, we need to store first ts in mem
        if not self._init_scope and self._frame == -1:
            self._frame += 1
            # can't simply return the same ts again- transformations would be applied twice
            # instead, return the pre-transformed copy
            return self._first_ts

        ts = self._read_frame(self._frame + 1)

        if self._init_scope:
            self._first_ts = self.ts.copy()
            self._init_scope = False

        return ts

    @property
    def n_frames(self):
        """Changes as stream is processed unlike other readers"""
        raise RuntimeError(
            "{}: n_frames is unknown".format(self.__class__.__name__)
        )

    def __len__(self):
        raise RuntimeError(
            "{} has unknown length".format(self.__class__.__name__)
        )

    def next(self):
        """Don't rewind after iteration. When _reopen() is called,
        an error will be raised
        """
        try:
            ts = self._read_next_timestep()
        except (EOFError, IOError):
            # Don't rewind here like we normally would
            raise StopIteration from None
        else:
            for auxname, reader in self._auxs.items():
                ts = self._auxs[auxname].update_ts(ts)

            ts = self._apply_transformations(ts)

        return ts

    def rewind(self):
        """Raise error on rewind"""
        raise RuntimeError(
            "{}: Stream-based readers can't be rewound".format(
                self.__class__.__name__
            )
        )

    # Incompatible methods
    def copy(self):
        raise NotImplementedError(
            "{} does not support copying".format(self.__class__.__name__)
        )

    def _reopen(self):
        if self._reopen_called:
            raise RuntimeError(
                "{}: Cannot reopen stream".format(self.__class__.__name__)
            )
        self._frame = -1
        self._reopen_called = True

    def __getitem__(self, frame):
        """Return the Timestep corresponding to *frame*.

        If `frame` is a integer then the corresponding frame is
        returned. Negative numbers are counted from the end.

        If frame is a :class:`slice` then an iterator is returned that
        allows iteration over that part of the trajectory.

        Note
        ----
        *frame* is a 0-based frame index.
        """
        if isinstance(frame, slice):
            _, _, step = self.check_slice_indices(
                frame.start, frame.stop, frame.step
            )
            if step is None:
                return FrameIteratorAll(self)
            else:
                return StreamFrameIteratorSliced(self, step)
        else:
            raise TypeError(
                "Streamed trajectories must be an indexed using a slice"
            )

    def check_slice_indices(self, start, stop, step):
        if start is not None:
            raise ValueError(
                "{}: Cannot expect a start index from a stream, 'start' must be None".format(
                    self.__class__.__name__
                )
            )
        if stop is not None:
            raise ValueError(
                "{}: Cannot expect a stop index from a stream, 'stop' must be None".format(
                    self.__class__.__name__
                )
            )
        if step is not None:
            if isinstance(step, numbers.Integral):
                if step < 1:
                    raise ValueError(
                        "{}: Cannot go backwards in a stream, 'step' must be > 0".format(
                            self.__class__.__name__
                        )
                    )
            else:
                raise ValueError(
                    "{}: 'step' must be an integer".format(
                        self.__class__.__name__
                    )
                )

        return start, stop, step

    def __getstate__(self):
        raise NotImplementedError(
            "{} does not support pickling".format(self.__class__.__name__)
        )

    def __setstate__(self, state: object):
        raise NotImplementedError(
            "{} does not support pickling".format(self.__class__.__name__)
        )

    def __repr__(self):
        return (
            "<{cls} {fname} with continuous stream of {natoms} atoms>"
            "".format(
                cls=self.__class__.__name__,
                fname=self.filename,
                natoms=self.n_atoms,
            )
        )


class StreamFrameIteratorSliced(FrameIteratorBase):

    def __init__(self, trajectory, step):
        super().__init__(trajectory)
        self._step = step

    def __iter__(self):
        # Calling reopen tells reader
        # it can't be reopened again
        self.trajectory._reopen()
        return self

    def __next__(self):
        try:
            # Burn the timesteps until we reach the desired step
            # Don't use next() to avoid unnecessary transformations
            while self.trajectory._frame + 1 % self.step != 0:
                self.trajectory._read_next_timestep()
        except (EOFError, IOError):
            # Don't rewind here like we normally would
            raise StopIteration from None

        return self.trajectory.next()

    def __len__(self):
        raise RuntimeError(
            "{} has unknown length".format(self.__class__.__name__)
        )

    def __getitem__(self, frame):
        raise RuntimeError("Sliced iterator does not support indexing")

    @property
    def step(self):
        return self._step
