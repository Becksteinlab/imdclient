Welcome to IMDClient's documentation!
=========================================================

IMDClient is a Python-based client library designed to interact with simulation engines making interactive molecular dynamics
(IMD) possible. It provides tools for communicating with the IMDv2 and new IMDv3 protocols, enabling seamless integration
between its API with simulation engines such as NAMD, LAMMPS, and GROMACS. Users can leverage this library to build
applications that interact with the IMD system efficiently, e.g., MDAnalysis IMDReader. IMDv3 has been standardized across
supported simulation engines with regard to transmitted data and is processed in IMDClient to provide consistent
`units <protocol_v3.html#units>`_. See `IMDv3 Protocol <protocol_v3.html>`_ for more details. IMDClient is ideal for
developers looking to produce on-the-fly analysis, avoiding storage of large trajectories.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   usage
   api
   protocol_v3
   header_def
   protocol_v2