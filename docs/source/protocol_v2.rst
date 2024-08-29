.. _v2_spec:

IMD Protocol in GROMACS
=======================

For developers interested in implementing the IMD v2 protocol, this document provides a fairly comprehensive specification of the protocol. 
The protocol is implemented in GROMACS in all versions starting from 5.0. The source code can be found
in the ``src/gromacs/imd`` directory.

- Version: 2
- Header size: 8 bytes (two 32 bit integers)
- Endianness: Big endian used for headers, except in "go" packet- see protocol steps below. Either endianness acceptable for packet bodies (position, force, and energy data).
- Number of clients per server: only 1 allowed

Header types
------------

Headers sent by the GROMACS server to the client
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 30 90
   :header-rows: 1

   * - Type
     - Enum Integer Value
     - Description
   * - IMD_ENERGIES
     - 1
     - Energy data for the system. For more information, see the `GROMACS code manual <https://manual.gromacs.org/5.1/doxygen/html-full/structIMDEnergyBlock.xhtml>`_
   * - IMD_FCOORDS
     - 2
     - Position data for all atoms in the system
   * - IMD_HANDSHAKE
     - 4
     - Sent to the client to inform the client of the version of the protocol being used as well as the endianness of the data

Headers sent by the client to the GROMACS server
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 30 90
   :header-rows: 1

   * - Type
     - Enum Integer Value
     - Description
   * - IMD_DISCONNECT
     - 0
     - Tells the GROMACS server to disconnect its client socket. If ``-imdwait`` is set, the GROMACS simulation
       will pause until the client reconnects and performs the handshake again.
   * - IMD_GO
     - 3
     - After receiving the handshake from GROMACS, the client has 1 second to send this signal to the server to begin
       receiving data. If the server does not receive this signal, it will disconnect the client and wait for another connection. If the GROMACS server
       was started with ``-imdwait``, the simulation will not begin until this signal is received.
   * - IMD_KILL
     - 5
     - Tells GROMACS to stop the simulation (requires that the simulation was started with ``-imdterm``)
   * - IMD_MDCOMM
     - 6
     - Inputs forces into the GROMACS simulation (requires that the simulation was started with ``-imdpull``)
   * - IMD_PAUSE
     - 7
     - After simulation has started, the client can send this signal to the GROMACS server to toggle the 
       simulation between paused and running states. This sends GROMACS into a blocking loop until
       the client sends another IMD_PAUSE signal.
   * - IMD_TRATE
     - 8
     - Adjusts the transmission rate in timesteps of the data from the server to the client. A value of 0 resets the rate to the default.

Unused headers
^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 30 90
   :header-rows: 1

   * - Type
     - Enum Integer Value
     - Description
   * - IMD_IOERROR
     - 9
     - Used internally by GROMACS to signal an error in the I/O operations with the client. This is not sent by the client.
   * - Count
     - 10
     - Defined in GROMACS source code but unused

Protocol steps
--------------

1. Pre-connection
^^^^^^^^^^^^^^^^^

.. code-block:: none

  GROMACS 1: Read -imdwait, -imdpull, and -imdterm options from the command line
  GROMACS 2: If -imdwait is set, block simulation run until GO signal received

2. Connection
^^^^^^^^^^^^^

.. code-block:: none

  GROMACS 1: Create a TCP socket and bind it to a port
  GROMACS 2: Listen for incoming connections, checking every 1 second
  Client 1: Connect to the server
  GROMACS 3: Accept the incoming connection, binding it to a new socket
  GROMACS 4: Send the handshake signal to the client:

    Header: 
        4 (int32) (IMD_HANDSHAKE)
        2 (int32) (Protocol version, must be unswapped so client can determine endianness)

  GROMACS 5: Wait up to 1 second to receive a GO signal from the client:

    Header:
        3 (int32) (IMD_GO)
        <val> (Unused attribute in header)

  * If the GO signal is not received, disconnect the client and wait for another connection *
  * After the client connects again, repeat the handshake and GO signal steps *

  Client 2: Send the GO signal to the server. All subsequent header packets will have big endianness 
            and all data packets will have the endianness specified in the handshake.
  
3. Simulation loop
^^^^^^^^^^^^^^^^^^

.. code-block:: none

  GROMACS (Server) 1: In every simulation integration step loop, check for incoming packets from the 
                      client. These packets can be one of:

    Header:
        5 (int32) (IMD_KILL)
        <val> (Unused attribute in header)

    Header:
        0 (int32) (IMD_DISCONNECT)
        <val> (Unused length attribute in header)

    Header:
        6 (int32) (IMD_MDCOMM)
        <val> (int32) (Number of forces that will be sent in the packet)

        Force packet:
            <val> (int32) (n indices of atoms to apply force to)
            <val> (float32) (n * 3 forces to apply to the atoms at the corresponding indices)

    Header:
        7 (int32) (IMD_PAUSE)
        <val> (Unused length attribute in header)

    Header:
        8 (int32) (IMD_TRATE)
        <val> (int32) (New transfer rate. Value of 0 means reset to default)

    * If any other header is recieved, disconnect the client and print an error message and wait for 
      another client connection *
    * After the client connects again, repeat the handshake and GO signal steps *

    GROMACS (Server) 2: Send energies and position data to the client if this timestep step lands
                        in the rate specified by the client (or the default, every timestep)

      Header:
          1 (int32) (IMD_ENERGIES)
          1 (int32) (Contstant value)

        Energy packet:
          <val> (float32) (1 float with the timestep and 9 floats describing the energy of the system)
      Header:
          2 (int32) (IMD_FCOORDS)
          <val> (int32) (Number of atoms in the system)

          Position packet:
            <val> (float32) (n atoms * 3 floats describing the position of each atom in the system)