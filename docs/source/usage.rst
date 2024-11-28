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
In GROMACS, you can use ``gmx mdrun`` with the ``-imdwait`` flag
to ensure that GROMACS will wait for a client before starting the simulation.
In GROMACS, you will know that the simulation is ready and waiting for the
MDAnalysis IMDReader client when this line is printed to the terminal:

.. code-block:: none

    IMD: Will wait until I have a connection and IMD_GO orders.

NAMD
----
In NAMD, the simulation will wait for a client connection when the  
``IMDon`` option is set to ``yes`` in the NAMD configuration file. 
Other options that can be set are detailed 
`here <https://github.com/amruthesht/namd-3.0/blob/IMDv3-dev/IMDv3-dev.md>`_. 
This will produce a simulation that is ready for a client connection with the 
following terminal message:

.. code-block:: none

    Info: INTERACTIVE MD AWAITING CONNECTION


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
