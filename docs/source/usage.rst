Streaming trajectories from a simulation engine with IMDv3
==========================================================

Configuring the simulation engine for IMDv3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To stream a trajectory from a simulation engine that supports IMDv3, 
first use the appropriate input options on the simulation engine 
to prepare it for the IMDClient receiver.

Below, we have proided brief instructions on how to setup the various 
simulation engine to output stream data using IMDv3.

GROMACS
-------
The IMDv3 protocol is currently not available as part of the official GROMACS release
or source code repository. However, the feature is currently available for use in 
the *imd-v3* branch of the forked repository
https://gitlab.com/heydenlabasu/streaming-md/gromacs/-/tree/imd-v3 ; clone the 
repository and build GROMACS from source.

In GROMACS, you can use ``gmx mdrun`` with the ``-imdwait`` flag
to ensure that GROMACS will wait for a client before starting the simulation.
Specific IMDv3 options were added to GROMACS, as documented in the
GROMACS User Guide `Molecular dynamics parameters (.mdp options)`_ under
*Interactive Molecular Dynamics (IMD)*.

In GROMACS, you will know that the simulation is ready and waiting for the
client when this line is printed to the terminal:

.. code-block:: none

    IMD: Will wait until I have a connection and IMD_GO orders.

You are now ready to connect to the simulation engine with a client.

.. TODO: update to official GROMACS docs (issue #79)
.. _`Molecular dynamics parameters (.mdp options)`:
   https://gitlab.com/heydenlabasu/streaming-md/gromacs/-/blob/imd-v3/docs/user-guide/mdp-options.rst?ref_type=heads&plain=1



NAMD
----
The IMDv3 protocol has been implementaed in NAMD and will be made available through the
official NAMD release in the near future. It is currently available as a part of the
official NAMD GitLab repository - https://gitlab.com/tcbgUIUC/namd.

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

You are now ready to connect to the simulation engine with a client.

LAMMPS
------
The IMDv3 protocol is part of the official LAMMPS distribution since ``patch_4Feb2025``.
It is available in the LAMMPS source code repository - https://github.com/lammps/lammps.
Information on using IMDv3 with LAMMPS can be found in the LAMMPS documentation
- https://docs.lammps.org/fix_imd.html

To use IMDv3 with LAMMPS, add the following lines to your LAMMPS input script:

.. code-block:: none

    fix ID group-ID imd <port> trate <frequency> version 3 unwrap <on/off> time <on/off> box <on/off> coordinates <on/off> velocities <on/off> forces <on/off>

Once the simulation is ready for a client connection, it will print 
following terminal message:

.. code-block:: none

    Waiting for IMD connection on port <port>

You are now ready to connect to the simulation engine with a client.

Using IMDClient with MDAnalysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the simulation is ready for a client connection, setup your MDAnalysis :class:`Universe`
like this: ::

    from IMDClient.IMD import IMDReader
    import MDAnalysis as mda
    # Pass host and port of the listening simulation
    # engine as the trajectory argument

    # GROMACS
    u = mda.Universe("topology.gro", "imd://localhost:8888")
    # NAMD
    u = mda.Universe("topology.psf", "imd://localhost:8888")

While this package allows the IMDReader to be automatically selected
based on the trajectory URL matching the pattern 'imd://<host>:<port>',
the format can be explicitly selected by passing the keyword argument
'format="IMD"' to the :class:`Universe`.
