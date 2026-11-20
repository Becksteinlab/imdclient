Streaming trajectories with IMDv2
=================================

IMDv2 is the long-standing interactive MD socket protocol used with VMD,
NAMD, GROMACS, and LAMMPS. IMDClient can connect to an IMDv2 simulation
and read coordinate (and optionally energy) streams over TCP.

For new projects, :doc:`usage` (IMDv3) is recommended: it has a published
specification, consistent wire units, and richer per-frame data (time,
box, velocities, forces). Use IMDv2 when your simulation is already set
up for the older protocol or your engine build only exposes v2.

Limitations of IMDv2
^^^^^^^^^^^^^^^^^^^^

There is no formal IMDv2 protocol document; behavior is inferred from
engine implementations. Relative to IMDv3, users should be aware of the
following:

* **Session metadata** — There is no handshake “session info” packet. The
  client cannot ask the engine which fields will appear in each frame.

* **Packet order** — Engines do not guarantee a fixed order of packet
  types per frame. Multiple energy packets can appear before coordinates
  (IMDClient keeps the latest energies and may warn).

* **Pause control** — IMDv2 uses a single pause toggle (``IMD_PAUSE``)
  instead of separate pause and resume packets.

* **After disconnect** — ``continue_after_disconnect`` and ``IMD_WAIT``
  apply only to IMDv3. Leave ``continue_after_disconnect`` unset
  (``None``) for v2 connections.

* **Transmission rate** — In LAMMPS (``trate``) and NAMD (``IMDfreq``),
  the frame interval is usually set in the simulation input. For **GROMACS
  IMDv2**, ``IMD-nst`` in the ``.mdp`` file is **not** supported (that
  option is for IMDv3). GROMACS v2 defaults to sending every integrator
  step until the client sends ``IMD_TRATE`` after the go packet; pass the
  desired interval as
  :paramref:`~imdclient.IMDClient.IMDClient.transmission_rate`.

Configuring the simulation for IMDv2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All **GROMACS, LAMMPS, or NAMD builds** support IMDv2 streaming.
IMDClient learns the protocol version from the handshake.

Typical usage input snippets for each engine are provided below.

**GROMACS** — Select an IMD group in the ``.mdp``, then run ``mdrun`` with
IMD enabled, for example::

    gmx mdrun -s topol.tpr -imdport 8888 -imdwait

**LAMMPS** — Use ``fix imd`` with the port and update interval, for
example::

    fix 2 all imd 8888 nowait off trate 10

**NAMD** — Enable IMD and set port, wait, and frequency in the config
file, for example::

    IMDon yes
    IMDport 8888
    IMDwait on
    IMDfreq 10

Once the simulation is waiting for a client, connect from your analysis
script or application.

Using IMDClient with IMDv2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The client-side Python workflow is as follows: open a connection, construct
:class:`~imdclient.IMDClient.IMDClient`, the client sends ``IMD_GO``, and
you read frames with :meth:`~imdclient.IMDClient.IMDClient.get_imdframe`.

.. note::

    **GROMACS (IMDv2-only builds)** — You cannot set the IMD frame interval
   with ``IMD-nst`` in the ``.mdp`` file. The GROMACS engine starts with 
   a default interval of **one integrator step** until the client sends an 
   ``IMD_TRATE`` packet explicitly.
   One can set the interval when creating the client via
   ``transmission_rate`` (integration steps between streamed frames). This
   sends a ``IMD_TRATE`` packet after the go handshake.

   .. code-block:: python

       client = IMDClient(
           host, port, n_atoms=n_atoms, transmission_rate=10
       )

   Omitting ``transmission_rate`` leaves the default of one step per
   frame, which can overload the socket and your analysis.

For **LAMMPS** and **NAMD** v2, the interval is usually set in the
simulation input (``trate`` or ``IMDfreq``). ``transmission_rate`` is
optional on the client in that case.

See :doc:`protocol_v3` for the IMDv3 specification and how it extends 
on IMDv2.