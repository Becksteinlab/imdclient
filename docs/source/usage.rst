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

As a part of the IMDClient package, we have provided a minimal implementation of a Reader 
Class viz. :class:`MinimalReader`, which can be used to read trajectory data from the socket
by leveraging :class:`IMDClient`. It however needs a function to process topology informaion 
to get the number of atoms in the simulation. We do so below by using the ``MDAnalysis`` 
library, which in principle can be done with any function that can read a topology file.

Once the simulation is ready for a client connection, setup your
:class:`MinimalReader` object like this: ::

    import MDAnalysis as mda
    from imdclient.tests.MinimalReader import MinimalReader
    # Pass host and port of the listening simulation
    # engine as an argument to the Reader

    # GROMACS
    u_mda = mda.Universe("topology.gro")
    n_atoms = u_mda.atoms.n_atoms
    u = MinimalReader("imd://localhost:8888", n_atoms=n_atoms)
    # NAMD
    u_mda = mda.Universe("topology.psf")
    n_atoms = u_mda.atoms.n_atoms
    u = MinimalReader("imd://localhost:8888", n_atoms=n_atoms)
    # LAMMPS
    u_mda = mda.Universe("topology.data")
    n_atoms = u_mda.atoms.n_atoms
    u = MinimalReader("imd://localhost:8888", n_atoms=n_atoms)

This example class can be used as a starting point to implement your own reader class to 
read trajectory data from the socket and generate on-the-fly simulation analysis.
