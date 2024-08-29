Proposed Protocol v3
====================

The suggested changes to the protocol are as follows:

1. Via an ``.mdp`` file setting, the user should be able to specify which of simulation box, positions, forces, and velocities are sent.
   These should be rate-adjustable like with xtc and trr output settings but should be set to every step by default.

2. The IMD_HANDSHAKE packet should have a 1-byte body which contains the configuration settings the simulation was setup with.
   This will allow the client to be choose appropriate packets to send and receive without redundant configuration.

   The modified handshake packet would look like this:

.. code-block:: none

    Header: 
        4 (int32) (IMD_HANDSHAKE)
        2 (int32) (Protocol version, must be unswapped so client can determine endianness)
    Body:
        <val> (bit) (imdpull: true or false. if true, the client can input forces to the simulation)
        <val> (bit) (imdwait: true or false. if true, the simulation will wait for the client to send a go signal before starting)
        <val> (bit) (imdterm: true or false. if true, the client can terminate the simulation)
        <val> (bit) (wrapped positions: true or false. if positions rate is 0, this is a placeholder value)

        <val> (uint32) (energies rate: number of steps that elapse between steps that contain energy data. 0 means never)
        <val> (uint32) (dimensions rate: number of steps that elapse between steps that contain dimension data. 0 means never)
        <val> (uint32) (positions rate: number of steps that elapse between steps that contain position data. 0 means never)
        <val> (uint32) (velocities rate: number of steps that elapse between steps that contain velocity data. 0 means never)
        <val> (uint32) (forces rate: number of steps that elapse between steps that contain force data. 0 means never)

   "wrapped positions" will be a new ``.mdp`` setting which specifies whether the atoms' positions
   should be adjusted to fit within the simulation box before sending. This is useful for visualization purposes.

3. The server should wait longer than 1 second (possibly up to 60s) for the go signal so that the client 
   has plenty of time to allocate memory buffers based on the endianness and information on included data types 
   it received in the handshake packet.

4. In the simulation loop, the server will send the client data in this order (if the configuration says to send it)
    
    i. Energy data (IMD_ENERGIES) unchanged
    
    ii. Dimension data (IMD_BOX) in triclinic vectors

    iii. Position data (IMD_FCOORDS) unchanged except box adjustments (see 5)
    
    iv. Velocity data (IMD_VELS) in the same manner as positions
    
    v. Force data (IMD_FORCES) in the same manner as positions

5. The server will send a new IMD_EOS (end of stream) packet after the last frame is sent unless the client initiates the disconnection with
   IMD_DISCONNECT.
