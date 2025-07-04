IMDv3 Protocol
==============

IMD (Interactive Molecular Dynamics) is a protocol for communicating molecular simulation 
data through a socket. 
It allows for two-way communication: via IMD, a simulation engine sends data to a receiver 
and the receiver can send forces and certain requests back to the simulation engine.

Version numbers in IMD are monotonically increasing integers. 
IMDv3, the protocol version described in this document, builds upon IMDv2, which is implemented 
at the time of writing in NAMD, GROMACS, and LAMMPS-4Feb2025. IMDv2 is itself 
a modification of the original IMD protocol published in 2001 :cite:p:`IMDv1:2001`.

.. list-table::
   :widths: 10 30
   :header-rows: 1

   * - IMD version
     - Protocol specification
   * - 1
     - *A system for interactive molecular dynamics simulation* :cite:p:`IMDv1:2001`
   * - 2
     - No official specification, but sample protocol and API implementation available :cite:p:`IMDv2`
   * - 3
     - This document

.. note:: 

   Because IMDv2 has no official specification, characteristics of the IMDv2 protocol mentioned in this document are inferred from 
   the sample implementation and commonalities in existing IMDv2 implementations.

Terms
-----
**Simulation engine**: The producer of simulation data which listens for a receiver connection and sends it simulation data.

**Receiver**: The receiver of simulation data, connects to listening simulation engine and receives simulation data through a socket.

**IMD frame**: A simulation integration step in which IMD data will be sent from the simulation engine to the receiver. The IMD 
transmission rate defines which integration steps are IMD frames.

**IMD transmission rate**: The number of integration steps in between each IMD frame. For example,
if the transmission rate is 1, every integration step in the simulation will be an IMD frame.

**IMD system**: The subset of all atoms in the simulation for which the simulation engine can output IMD data
and for which the receiver can send forces to the simulation engine to apply to.

**IMD energy block**: A data structure specific to IMD which contains the simulation integration step and information
on the energy of the simulated system (not just the IMD system).

Example
-------

For a quick understanding of how IMDv3 works, this scenario describes the flow of a typical IMDv3 session:

1. GROMACS, a simulation engine, starts listening for IMDv3 connections on port 8888
2. IMDClient, an IMDv3 receiver, connects to the listening port
3. Through the connected socket, GROMACS sends IMDClient a handshake packet containing its running machine's endianness and the version of the session (v3)
4. GROMACS sends IMDClient a session info packet (introduced in IMDv3) informing it that it should expect time, box, coordinate, and velocity data in this session
5. After parsing the handshake and session info packets, IMDClient sends a go signal to GROMACS and begins waiting for data packets
6. GROMACS repeatedly sends time, box, coordinate, and velocity packets for each simulation frame in a fixed order
7. IMDClient continuously parses the stream of bytes, making the data available for processing by analysis code
8. GROMACS receives a termination signal from the environment it is running in and closes the socket
9. IMDClient recognizes that the session has ended, raising an EOF error after the last available frame of data has been processed

Packet Types
------------

In IMD, there are two kinds of packets: header packets and body packets.

A header packet is composed of 8 bytes. The first 4 bytes 
are always the header type and the next 4 bytes 
serve as a flexible slot for holding other information.
All bytes in header packets are provided in network (big)
endian order except for the second 4 
bytes of the handshake packet. This is described in greater detail
in the :ref:`handshake <handshake>` section.

In this document, header packets are described like this:

.. code-block:: none

   Header:
      <val> (int32) Header type
      <val> (dtype) Data slot description

Some header packets sent by the simulation engine have associated body packets that contain additional data, 
like atomic coordinates. These body packets vary in length and composition,
but their bytes are always encoded in the native endianness of the machine running the simulation engine.

This table lists all header types available in IMDv3. Below, each header type 
and its associated body packet (if present) is described in detail.

.. list-table::
   :widths: 10 10 10 10
   :header-rows: 1

   * - Header type
     - 32-bit integer enum value
     - In IMDv2
     - In IMDv3
   * - :ref:`disconnect`
     - 0
     - ✅
     - ✅
   * - :ref:`energies`
     - 1
     - ✅
     - ✅
   * - :ref:`coordinates`
     - 2
     - ✅
     - ✅
   * - :ref:`go`
     - 3 
     - ✅
     - ✅
   * - :ref:`handshake`
     - 4
     - ✅
     - ✅
   * - :ref:`kill`
     - 5
     - ✅
     - ✅
   * - :ref:`md-communication`
     - 6
     - ✅
     - ✅
   * - :ref:`pause`
     - 7 
     - ✅
     - ✅
   * - :ref:`transmission-rate`
     - 8
     - ✅
     - ✅
   * - :ref:`io-error`
     - 9
     - ✅
     - ✅
   * - :ref:`session-info`
     - 10
     - ❌
     - ✅
   * - :ref:`resume`
     - 11
     - ❌
     - ✅
   * - :ref:`time`
     - 12
     - ❌
     - ✅
   * - :ref:`box`
     - 13
     - ❌
     - ✅
   * - :ref:`velocities`
     - 14
     - ❌
     - ✅
   * - :ref:`forces`
     - 15
     - ❌
     - ✅
   * - :ref:`wait`
     - 16
     - ❌
     - ✅

Packet Header Definitions
-------------------------

.. _disconnect:

Disconnect
^^^^^^^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to indicate that the simulation engine should 
close the connected socket. Whether the simulation engine pauses execution until another connection is
made is an implementation decision.

.. code-block:: none

   Header:
      0 (int32) Disconnect
      <val> (no type) Unused slot, any value acceptable

.. _energies:

Energies
^^^^^^^^

Sent from the simulation engine to the receiver each IMD frame if 
energies were previously specified for this session in the :ref:`session info packet <session-info>`.

.. note:: 
  While the integration step is included in this
  packet, this is a result of inheriting the IMD energy block from IMDv2. It is recommended
  to make use of the 64-bit integer integration step value from the :ref:`time packet <time>`
  in analysis code instead.

.. code-block:: none

   Header:
      1 (int32) Energies
      1 (int32) Number of IMD energy blocks being sent

   Body:
      <val> (int32) Current integration step of the simulation
      <val> (float32) Absolute temperature
      <val> (float32) Total energy
      <val> (float32) Potential energy
      <val> (float32) Van der Waals energy
      <val> (float32) Coulomb interaction energy
      <val> (float32) Bonds energy
      <val> (float32) Angles energy
      <val> (float32) Dihedrals energy
      <val> (float32) Improper dihedrals energy

.. _coordinates:

Coordinates
^^^^^^^^^^^

Sent from the simulation engine to the receiver each IMD frame if 
coordinates were previously specified for this session in the :ref:`session info packet <session-info>`.

.. code-block:: none

   Header:
      2 (int32) Coordinates
      <n_atoms> (int32) Number of atoms in the IMD system

   Body:
      <array> (float32[n_atoms * 3]) X, Y, and Z coordinates of each atom in the 
                                     IMD system encoded in the order 
                                     [X1, Y1, Z1, ..., Xn, Yn, Zn]

.. _go:

Go
^^

Sent from the receiver to the simulation engine after the receiver receives 
the :ref:`handshake <handshake>` and :ref:`session info <session-info>` packets. 

If the simulation engine does not 
receive this packet within 1 second of sending the handshake and session info 
packets, it should assume the receiver is incompatible. Whether the simulation engine
exits or accepts another connection after this is an implementation decision.

.. code-block:: none

   Header:
      3 (int32) Go
      <val> (no type) Unused slot, any value acceptable

.. _handshake:

Handshake
^^^^^^^^^

Sent from the simulation engine to the receiver after a socket connection
is established. Unlike other header packets, the last four bytes of this packet are provided in 
the native endianness of the sending simulation engine's hardware.

The receiver can use this packet to determine both the IMD version
of the session and the endianness of the simulation engine. By providing 
the endianness of the machine running the simulation engine, the bulk of the 
data being sent in the session, i.e. the body packets, do not have to be swapped 
by the simulation engine before being sent, speeding up execution.

.. code-block:: none

   Header:
      4 (int32) Handshake
      3 (int32, unswapped byte order) IMD version used in session

.. _kill:

Kill
^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to request that the simulation engine
stops execution of the simulation and exits. Whether or not the simulation engine 
honors this request is an implementation decision.

.. code-block:: none

   Header:
      5 (int32) Kill
      <val> (no type) Unused slot, any value acceptable

.. _md-communication:

MD Communication
^^^^^^^^^^^^^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to request that the forces 
in the body packet are applied to the atoms specified in the body packet. 
Whether or not the simulation engine honors this request is an implementation decision.

.. code-block:: none

   Header:
      6 (int32) MD Communication
      <n_atoms> (int32) Number of atoms in the IMD system to apply forces to

   Body:
      <array> (int32[n_atoms]) Indices of atoms in the IMD system to apply forces to
      <array> (float32[n_atoms * 3]) The X, Y, and Z components of forces to be applied to
                                     the atoms at the indices specified in the above array
     
The array of IMD system indices does not need to be monotonically increasing, meaning 
the indices can be "out of order". However, the index array cannot contain any index twice. 
Force vectors acting on the same index should 
be combined before being sent to the simulation engine to be applied.

.. note:: 
   
   Though this packet is sent by the receiver, the rule that all body packets are 
   sent in the native endianness of the machine running the simulation engine
   still applies here. The receiver must use the endianness it gets from 
   the :ref:`handshake <handshake>` and swap the endianness of the indices and forces 
   if necessary before sending.

.. _pause:

Pause
^^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to request that the simulation
engine pauses execution of the simulation until a :ref:`resume packet <resume>` is sent.
Pause is idempotent, meaning subsequent pause packets sent after the first one will have no effect.


.. code-block:: none

   Header:
      7 (int32) Pause
      <val> (no type) Unused slot, any value acceptable

.. versionchanged:: 3

   In IMDv2, pause acted as a toggle, meaning sending a pause packet twice 
   would pause and then resume the simulation's execution. In IMDv3, the :ref:`resume packet <resume>`
   is required to resume a paused simulation since pausing is idempotent.

.. _transmission-rate:

Transmission rate
^^^^^^^^^^^^^^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to change the IMD transmission rate. 

.. code-block:: none

   Header:
      8 (int32) Transmission rate
      <val> (int32) New transmission rate. Any value less than 1 will reset 
                    the transmission rate to its default value (configured
                    by the simulation engine)

.. _io-error:

IO Error
^^^^^^^^

Never sent from one party to another during an IMD session. Can be used internally 
by the simulation engine or receiver to indicate an error has occurred.

.. code-block:: none

   Header:
      9 (int32) IO Error
      <val> (no type) Unused slot, any value acceptable

.. _session-info:

Session info
^^^^^^^^^^^^

Sent by the simulation engine to the receiver immediately after
the :ref:`handshake` is sent to indicate to the receiver which data it 
should expect for each IMD frame during the session along with
whether coordinates will be wrapped into the simulation box if present.

.. code-block:: none

   Header:
      10 (int32) Session info
      7 (int32) Number of 1-byte configuration options in the body packet
    
   Body:
      <val> (int8) Nonzero if time packets sent in each IMD frame
      <val> (int8) Nonzero if IMD energy block packets sent in each IMD frame
      <val> (int8) Nonzero if box packets sent in each IMD frame
      <val> (int8) Nonzero if coordinate packets sent in each IMD frame 
      <val> (int8) Nonzero if coordinates wrapped into the simulation box. 
                   Meaningless if coordinates not sent in the session
      <val> (int8) Nonzero if velocity packets sent in each IMD frame 
      <val> (int8) Nonzero if force packets sent in each IMD frame 

.. versionadded:: 3

.. _resume:

Resume
^^^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to request that the simulation resumes execution
if it is in a paused state. Like :ref:`pause <pause>`, resume is idempotent.

.. code-block:: none

   Header:
      11 (int32) Resume
      <val> (no type) Unused slot, any value acceptable

.. versionadded:: 3

.. _time:

Time
^^^^

Sent from the simulation engine to the receiver each IMD frame if 
time packets were previously specified for this session in the :ref:`session info packet <session-info>`.

.. code-block:: none

   Header:
      12 (int32) Time
      1 (int32) Number of time packets being sent

   Body:
      <val> (float64) Simulation integrator time step (dt)
      <val> (float64) Current simulation time
      <val> (int64) Current integration step of the simulation

.. versionadded:: 3

.. _box:

Box
^^^

Sent from the simulation engine to the receiver each IMD frame if 
box packets were previously specified for this session in the :ref:`session info packet <session-info>`.

.. code-block:: none

   Header:
      13 (int32) Box
      1 (int32) Number of simulation boxes being sent
   Body:
      <array> (float32[9]) Triclinic box vectors for the simulation encoded in 
                           in the order [ABC] where A = (aX,aY,aZ), B = (bX,bY,bZ), 
                           and C = (cX,cY,cZ)

.. versionadded:: 3

.. _velocities:

Velocities
^^^^^^^^^^

Sent from the simulation engine to the receiver each IMD frame if 
velocities were previously specified for this session in the :ref:`session info packet <session-info>`.

.. code-block:: none

   Header:
      14 (int32) Velocities
      <n_atoms> (int32) Number of atoms in the IMD system

   Body:
      <array> (float32[n_atoms * 3]) X, Y, and Z components of the velocities 
                                     of each atom in the 
                                     IMD system encoded in the order 
                                     [Vx1, Vy1, Vz1, ..., Vxn, Vyn, Vzn]

.. versionadded:: 3

.. _forces:

Forces
^^^^^^

Sent from the simulation engine to the receiver each IMD frame if 
forces were previously specified for this session in the :ref:`session info packet <session-info>`.

.. code-block:: none

   Header:
      15 (int32) Forces
      <n_atoms> (int32) Number of atoms in the IMD system

   Body:
      <array> (float32[n_atoms * 3]) X, Y, and Z components of the forces 
                                     of each atom in the 
                                     IMD system encoded in the order 
                                     [Fx1, Fy1, Fz1, ..., Fxn, Fyn, Fzn]

.. versionadded:: 3

.. _wait:

Wait
^^^^

Sent from the receiver to the simulation engine any time after the :ref:`session info packet <session-info>`
has been sent to request that the simulation engine modify its waiting behavior mid-simulation either
from blocking to non-blocking or vice versa.
Whether or not the simulation engine honors this request is an implementation decision. 

Regardless of whether this packet is accepted, the simulation engine will have an initial waiting behavior which applies
to the beginning of the simulation:

1. Blocking: Wait until a receiver is connected to begin execution of the simulation 
2. Non-blocking: Begin the simulation regardless of whether a receiver is connected and continuously check on the listening socket for a receiver attempting to connect 

The simulation engine's waiting behavior also applies when a receiver disconnects mid-simulation:

1. Blocking: Pause simulation execution and wait until a receiver is connected to resume execution 
2. Non-blocking: Continue execution, continuously checking on the listening socket for a receiver attempting to connect

 .. code-block:: none

   Header:
      16 (int32) Wait
      <val> (int32) Nonzero to set the simulation engine's waiting behavior to blocking, 0
                    to set the simulation engine's waiting behavior to non-blocking

.. note:: 

   The purpose of this packet is to allow a receiver to monitor the first *n* frames 
   of a simulation and then disconnect without blocking the continued execution of the 
   simulation.

.. versionadded:: 3

Packet order
------------

After the simulation engine sends the :ref:`handshake <handshake>` and :ref:`session info <session-info>`
packets to the receiver and gets back a :ref:`go <go>` signal, it begins sending simulation data via
IMD. The data within each IMD frame is always sent in the same, fixed order:

1. Time
2. Energy block
3. Box
4. Coordinates
5. Velocities
6. Forces

If the simulation engine is configured to send only a strict subset of all
available data packets, the fixed order of the list still applies to the
remaining packets in the session. 

.. versionchanged:: 3

   In IMDv2, any packet order sent by the simulation engine is acceptable
   and IMD frames in the same session don't have to contain the same data packets.
   For example, an IMD frame in which only energies are sent can be followed by 
   an IMD frame in which only coordinates are sent. In IMDv3, all packets 
   specified in the session info must be sent for every IMD frame and in the same order.

Units
-----

The units in IMDv3 are fixed. The simulation engine must convert 
values into these units before sending them through the socket. 
The receiver must also convert forces it sends back to the simulation 
engine into these units.


.. list-table::
   :widths: 10 10
   :header-rows: 1

   * - Measurement
     - Unit
   * - Length
     - angstrom
   * - Velocity
     - angstrom/picosecond
   * - Force
     - kilojoules/(mol*angstrom)
   * - Time
     - picosecond
   * - Energy
     - kilojoules/mol

IMD port number
---------------

The preferred port for IMD communication is 8888, but the simulation engine may freely specify the 
port at which it listens for a receiver.

References
----------

.. bibliography::
