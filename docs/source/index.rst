imdclient documentation
=======================

:Release: |release|
:Version: |version|
:Date: |today|
:License: |license|

:mod:`imdclient` is a Python-based client library designed to interact with
simulation engines that use the *interactive molecular dynamics* (IMD) protocol
for **streaming trajectory data while a simulation is running**.

The package supports the legacy IMDv2 protocol and the new IMDv3 protocol,
enabling seamless integration between its API with simulation engines such as
NAMD_, LAMMPS_, and GROMACS_. Users can leverage this library to build
applications that interact with the IMD system efficiently, e.g., MDAnalysis_
uses *imdclient* to work with IMD streams.

IMDv3 has been standardized across supported simulation engines with regard to
transmitted data and is processed in *imdclient* to provide consistent `units
<protocol_v3.html#units>`_; see `IMDv3 Protocol <protocol_v3.html>`_ for more
details. *imdclient* is ideal for developers looking to produce on-the-fly
analysis, avoiding storage of large trajectories.

*imdclient* also supports the legacy protocol version IMDv2; setup, usage, and
limitations of the protocol are described under :doc:`usage_v2`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   usage
   usage_v2   
   api
   protocol_v3


.. _NAMD: http://www.ks.uiuc.edu/Research/namd
.. _LAMMPS: https://www.lammps.org/
.. _GROMACS: https://www.gromacs.org
.. _MDAnalysis: https://www.mdanalysis.org
