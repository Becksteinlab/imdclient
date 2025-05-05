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

.. _Packet Header Definitions: header_def.html

See `Packet Header Definitions`_ for more information.

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
