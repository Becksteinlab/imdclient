imdclient documentation
=======================

:Release: |release|
:Version: |version|
:Date: |today|
:License: |license|
:Source code: `github.com/Becksteinlab/imdclient`_

:mod:`imdclient` is a Python-based client library designed to interact with
simulation engines that use the *interactive molecular dynamics* (IMD) protocol
for **streaming trajectory data while a simulation is running**.  It enables
on-the-fly analysis of MD simulations at much higher frequency than feasible
with typical post-simulation workflows as it avoids the high storage demands of
trajectories written at high output frequency. Other applications include
monitoring and visualization of running simulations, writing of arbitrary
simulation subsystems to trajectories at high frequency, workflows that
exchange trajectory data without requiring slow disk storage, and online
learning workflows.

The package supports the legacy IMDv2 protocol and the new IMDv3 protocol,
enabling seamless integration between its API with simulation engines such as
NAMD_, LAMMPS_, and GROMACS_. Users can leverage this library to build
applications that interact with the IMD system efficiently, e.g., MDAnalysis_
uses *imdclient* to work directly with IMD streams.

**IMDv3** has been standardized across supported simulation engines with regard to
transmitted data and is processed in *imdclient* to provide consistent `units
<protocol_v3.html#units>`_; see `IMDv3 Protocol <protocol_v3.html>`_ for more
details. *imdclient* also supports the legacy protocol version **IMDv2** for MD
engines that have not yet incorporated changes for IMDv3, even though the
protocol imposes severe limitations on trajectory analysis.

The package was primarily designed for analysis and not for the "interactive"
component of IMD, and thus injecting forces back into the MD engine is
currently not supported (see `issue GH-45`_).

.. SeeAlso::
   * `mdadash`_ is a browser-based real-time dashboard for monitoring, tracking
     and analyzing running MD simulations, based on MDAnalysis_ and *imdclient*.
   * MDAnalysis `IMDReader`_ uses *imdclient* to make IMD streams readable like
     standard trajectories (with some inherent limitations due to streaming)


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   usage
   usage_v2   
   api
   protocol_v3

.. _`github.com/Becksteinlab/imdclient`:
   https://github.com/Becksteinlab/imdclient
.. _NAMD: http://www.ks.uiuc.edu/Research/namd
.. _LAMMPS: https://www.lammps.org/
.. _GROMACS: https://www.gromacs.org
.. _MDAnalysis: https://www.mdanalysis.org
.. _`issue GH-45`: https://github.com/Becksteinlab/imdclient/issues/45 
.. _mdadash: https://mdadash.readthedocs.io
.. _IMDReader: https://docs.mdanalysis.org/stable/documentation_pages/coordinates/IMD.html
