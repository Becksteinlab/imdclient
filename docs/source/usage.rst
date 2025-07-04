Streaming trajectories from a simulation engine with IMDv3
==========================================================

Configuring the simulation engine for IMDv3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To stream a trajectory from a simulation engine that supports IMDv3, 
first use the appropriate input options on the simulation engine 
to prepare it for the IMDClient receiver.

Below, we have provided brief instructions on how to setup the various 
simulation engine to output stream data using IMDv3.

GROMACS
-------
In GROMACS, you can use ``gmx mdrun`` with the ``-imdwait`` flag
to ensure that GROMACS will wait for a client before starting the simulation.
In GROMACS, you will know that the simulation is ready and waiting for the
IMDClient when this line is printed to the terminal:

.. code-block:: none

    IMD: Will wait until I have a connection and IMD_GO orders.

You are now ready to connect to the simulation engine with a client.

NAMD
----
To use IMDv3 with NAMD, add the following lines to your NAMD configuration file:

.. code-block:: none

    IMDon yes
    IMDport <port, must be the same port used for the client>
    IMDwait <yes/no>
    IMDfreq <frequency of sending data to the client>

    IMDsendPositions <yes/no>
    IMDsendEnergies <yes/no>
    IMDsendTime <yes/no>
    IMDsendBoxDimensions <yes/no>
    IMDsendVelocities <yes/no>
    IMDsendForces <yes/no>
    IMDwrapPositions <yes/no>


Once the simulation is ready for a client connection, it will print 
following terminal message:

.. code-block:: none

    Info: INTERACTIVE MD AWAITING CONNECTION

You are now ready to connect to the simulation engine with the IMDClient.

LAMMPS
------
To use IMDv3 with LAMMPS, add the following lines to your LAMMPS input script:

.. code-block:: none

    fix ID group-ID imd <port> trate <frequency> version 3 unwrap <on/off> time <on/off> box <on/off> coordinates <on/off> velocities <on/off> forces <on/off>

Once the simulation is ready for a client connection, it will print 
following terminal message:

.. code-block:: none

    Waiting for IMD connection on port <port>

You are now ready to connect to the simulation engine with a client.

Using IMDClient
^^^^^^^^^^^^^^^

Once the simulation is ready for a client connection, one can setup
the client using the :class:`~imdclient.IMDClient` class: ::

    from imdclient.utils import parse_host_port
    from imdclient.IMDClient import IMDClient

    host, port = parse_host_port("imd://localhost:8888")

    # `n_atoms` is the number of atoms in the simulation
    # Adjust this value according to your simulation setup

    # This forms the connection and starts the simulation 
    # by sending the `IMD_GO`
    client = IMDClient(host, port, n_atoms=1000)

    # Read trajectory data from the IMDBuffer which stores
    # data received from the socket

    i = 0
    while True:
        try:
            frame = client.get_imdframe()
        except EOFError:
            break
        else:
            i += 1
            # Process and analyze the frame data as needed
            # For example, print frame number, simulation time, and positions of atom 0
            print(f"Frame {i}: time={frame.time}, atom 0 position={frame.positions[0]}")

The :meth:`~imdclient.IMDClient.get_imdframe` method returns an 
:class:`~imdclient.IMDFrame` object containing the frame data 
read from the buffer and received from the socket.

The above example can be used as a starting point to implement your own reader 
class that utilizes :class:`~imdclient.IMDClient` to read trajectory data 
from the socket and generate on-the-fly simulation analysis.

.. SeeAlso::
    `MDAnalysis <https://www.mdanalysis.org>`_ (from release 2.10.0 onwards) can
    directly read IMDv3 streams.
