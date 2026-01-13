.. _quickstart:

Quick Start Guide
=================

This guide will get you up and running with sandlersteam in just a few minutes.

Basic Concepts
--------------

sandlersteam provides implementations of steam-table lookups
from tables in Appendix III of Sandler's textbook.  These tables
report temperatures in °C and pressures in either kPa or MPa.  Queries
can be made using any two independent state variables, such as temperature
and pressure, or temperature and specific volume.  The following properties
are returned:

* Specific enthalpy (h), entropy (s), volume (v), and internal energy (u) 
* Vapor fraction (x\ :sub:`vap`) for two-phase states

Temperature is in °C and pressure is in megapascals (MPa) by default.  
Other units of pressure are supported.  
See the :doc:`API` documentation for details.

First Calculation
-----------------

Let's calculate the molar volume of methane at 400 K and 0.5 MPa using the Peng-Robinson equation:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

Unsaturated steam/water state calculations:

.. code-block:: bash

   sandlersteam state -T 400 -P 1.1

Output::

   T  =   400 C
   P  =   1.1 MPa
   v  =  0.2807 m3/kg
   s  =  7.42125 kJ/kg-K
   h  =  3262.3 kJ/kg
   u  =  2956.1 kJ/kg
   Pv =  308.77 kJ/kg


Saturated steam/water state calculations:

.. code-block:: bash

   sandlersteam state -P 4.0 -x 0.9

Output::

   T   =  250.4 C
   P   =     4 MPa
   v   =  0.0449272 m3/kg
   s   =  5.74273 kJ/kg-K
   h   =  2629.99 kJ/kg
   u   =  2450.3 kJ/kg
   Pv  =  179.709 kJ/kg
   x   = 0.9 kg vapor/kg total
   vL  =  0.001252 m3/kg
   sL  =  2.7964 kJ/kg-K
   hL  =  1087.31 kJ/kg
   uL  =  1082.31 kJ/kg
   PvL =  5.008 kJ/kg
   vV  =  0.04978 m3/kg
   sV  =  6.0701 kJ/kg-K
   hV  =  2801.4 kJ/kg
   uV  =  2602.3 kJ/kg
   PvV =  199.12 kJ/kg

From Python
~~~~~~~~~~~

.. code-block:: python

   from sandlersteam.state import State
   state1 = State(T=200.0, P=0.1).lookup()

   print(f"Specific volume: {state1.v} {state1.volume_units}/{state1.mass_units}")   


State-Change Calculations
-------------------------

Calculate property changes between two thermodynamic states:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sandlersteam delta -P1 4.0 -T1 800 -P2 4.5 -T2 1200 --show-states

Output::

   State-change calculations for water/steam:

   State 1:                       State 2:
   T  =   800 C                   T  =  1200 C
   P  =     4 MPa                 P  =   4.5 MPa
   v  =  0.12287 m3/kg            v  =  0.15098 m3/kg
   s  =  7.8502 kJ/kg-K           s  =  8.5825 kJ/kg-K
   h  =  4141.5 kJ/kg             h  =  5136.9 kJ/kg
   u  =   3650 kJ/kg              u  =  4457.5 kJ/kg
   Pv =  491.48 kJ/kg             Pv =  679.41 kJ/kg

   Property changes:
   ΔT  =   400 C
   ΔP  =   0.5 MPa
   Δv  =  0.02811 m3/kg
   Δs  =  0.7323 kJ/kg-K
   Δh  =  995.4 kJ/kg
   Δu  =  807.5 kJ/kg
   ΔPv =  187.93 kJ/kg

From Python
~~~~~~~~~~~

State-change calculations can be performed by creating two EOS objects representing the initial and final states, then using the built-in methods to compute property differences:

.. code-block:: python

   from sandlersteam.state import State
   
   # State 1
   state1 = State(T=350, P=7.5).lookup()
   
   # State 2
   state2 = State(T=400, P=15.5).lookup()
   
   # Calculate property differences using built-in methods
   deltas = state1.delta(state2)
   
   print(f"Δh = {deltas['h']:7g} {state1.energy_unit}/{state1.mass_unit}")
   print(f"Δs = {deltas['s']:7g} {state1.energy_unit}/{state1.mass_unit}-K")
   print(f"Δu = {deltas['u']:7g} {state1.energy_unit}/{state1.mass_unit}")

Next Steps
----------

* Learn more about the :doc:`cli` for advanced command-line usage
* Explore :doc:`examples` for more complex calculations
* Check the :doc:`api/API` for complete API documentation
