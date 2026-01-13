.. _examples:

Examples
========

This page provides practical examples of using sandlersteam for common thermodynamic calculations.

Basic Property Calculations
----------------------------

Example 1: Single State Point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate properties of superheated steam at 350 C and 2.0 MPa:

.. code-block:: python

   from sandlersteam import State
   
   # Initialize EOS and set state
   state = State(T=350, P=2.0).lookup()
   
   # Display results
   print(f"State conditions:")
   print(f"  T = {state.T} {state.temperature_units}")
   print(f"  P = {state.P} {state.pressure_units}")
   print(f"\nCalculated properties:")
   print(f"  v = {state.v:.6f} {state.volume_units}/{state.mole_units}")
   print(f"  h = {state.h:.2f} {state.energy_units}/{state.mole_units}")
   print(f"  s = {state.s:.4f} {state.energy_units}/{state.mole_units}-K")

Example 2: Saturated State
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate properties of saturated steam at 0.5 MPa:

.. code-block:: python

   from sandlersteam import State
   
   # Initialize EOS and set state
   state = State(P=0.5, x=1.0).lookup()  # x=1.0 for saturated vapor
   
   # Display results
   print(f"Saturated steam at P = {state.P} {state.pressure_units}:")
   print(f"  T = {state.T} {state.temperature_units}")
   print(f"  v = {state.v:.6f} {state.volume_units}/{state.mole_units}")
   print(f"  h = {state.h:.2f} {state.energy_units}/{state.mole_units}")
   print(f"  s = {state.s:.4f} {state.energy_units}/{state.mole_units}-K")

Property Changes
----------------

Example 3: Isothermal Compression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate property changes during isothermal compression:

.. code-block:: python

    import numpy as np
    from sandlersteam import State
   
    state1 = State(T=300, P=1.0).lookup()
    state2 = State(T=300, P=10.0).lookup()
    deltas = state1.delta(state2)

    print(f"Isothermal Compression of Superheated Steam at T = {state1.T} {state1.temperature_units}")
    print(f"P: {state1.P} to {state2.P} {state1.pressure_units}\n")
    print(f"State Changes:")
    print(f"  v: {state1.v:.6f} → {state2.v:.6f} {state1.volume_units}/{state1.mass_units}\n")
    print(f"Property Changes:")
    print(f"  Δh  = {deltas['h']:.2f} {state1.energy_units}/{state1.mass_units}")
    print(f"  Δs  = {deltas['s']:.4f} {state1.energy_units}/{state1.mass_units}-K")
    print(f"  Δu  = {deltas['u']:.4f} {state1.energy_units}/{state1.mass_units}")
    print(f"  ΔPv = {deltas['Pv']:.4f} {state1.energy_units}/{state1.mass_units}")

If we assume reversibility, we can also calculate the work done during compression:

.. code-block:: python

    Q = state1.T * deltas['s']  # Heat transfer
    W = deltas['h'] - Q        # Work done
    print(f"\nWork done during isothermal compression: {W:.2f} {state1.energy_units}/{state1.mass_units}")

Example 4: Isentropic Expansion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate work and property changes during isentropic expansion:

.. code-block:: python

    from sandlersteam import State
   
    state1 = State(T=500, P=10.0).lookup()
    s_target = state1.s  # Isentropic process
    state2 = State(s=s_target, P=1.0).lookup()
    deltas = state1.delta(state2)

    print(f"Isentropic Expansion of Superheated Steam from P = {state1.P} to {state2.P} {state1.pressure_units}\n")
    print(f"State Changes:")
    print(f"  T: {state1.T} → {state2.T} {state1.temperature_units}")
    print(f"  v: {state1.v:.6f} → {state2.v:.6f} {state1.volume_units}/{state1.mass_units}\n")
    print(f"Property Changes:")
    print(f"  Δh  = {deltas['h']:.2f} {state1.energy_units}/{state1.mass_units}")
    print(f"  Δs  = {deltas['s']:.4f} {state1.energy_units}/{state1.mass_units}-K (should be 0 for isentropic)")
    print(f"  Δu  = {deltas['u']:.4f} {state1.energy_units}/{state1.mass_units}")
    print(f"  ΔPv = {deltas['Pv']:.4f} {state1.energy_units}/{state1.mass_units}")

    W = deltas['h']  # Work done by the system
    print(f"\nWork delivered during isentropic expansion: {W:.2f} {state1.energy_units}/{state1.mass_units}")
