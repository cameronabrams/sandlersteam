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

Note that sandlersteam represents all dimensioned properties using :class:`pint.Quantity` objects, so you can easily convert between different units in the API.


First Calculation
-----------------

Let's calculate the molar volume of methane at 400 K and 0.5 MPa using the Peng-Robinson equation:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

Unsaturated steam/water state calculations:

.. code-block:: bash

   sandlersteam state -T 400 -P 1.1

Output::

   T  =   400 degC
   P  =   1.1 MPa
   v  =  0.2807 m**3 / kg
   s  =  7.42125 kJ / (kg * K)
   h  =  3262.3 kJ / kg
   u  =  2956.1 kJ / kg
   Pv =  308.77 kJ / kg


Saturated steam/water state calculations:

.. code-block:: bash

   sandlersteam state -P 4.0 -x 0.9

Output::

   T   =  250.4 degC
   P   =     4 MPa
   v   =  0.0449272 m**3 / kg
   s   =  5.74273 kJ / (kg * K)
   h   =  2629.99 kJ / kg
   u   =  2450.3 kJ / kg
   Pv  =  179.709 kJ / kg
   x   = 0.9 mass fraction vapor
   vL  =  0.001252 m**3 / kg
   sL  =  2.7964 kJ / (kg * K)
   hL  =  1087.31 kJ / kg
   uL  =  1082.31 kJ / kg
   PvL =  5.008 kJ / kg
   vV  =  0.04978 m**3 / kg
   sV  =  6.0701 kJ / (kg * K)
   hV  =  2801.4 kJ / kg
   uV  =  2602.3 kJ / kg
   PvV =  199.12 kJ / kg

From Python
~~~~~~~~~~~

.. code-block:: python

   from sandlersteam.state import State
   state1 = State(T=200.0, P=0.1)

   print(f"Specific volume: {state1.v}")

When using the API, note that you can use dot notation to access properties directly, and to reset the initial inputs. For example:

.. code-block:: python

   from sandlersteam.state import State
   state = State(T=200.0, P=0.1)

   print(f"Specific volume: {state.v}")

   # Change state to new conditions
   state.T = 300.0  # this will fully specify the state and recalculate all properties, at T=300 C and previous P=0.1 MPa
   print(f"New specific volume: {state.v}")

   state.P = 0.5 # this will fully specify the state and recalculate all properties at P=0.5 MPa and previous T=300 C
   print(f"Another new specific volume: {state.v}")

If you use dot notation to overwrite a property that is not a state's input property, this will wipe out the existing inputs and you will need to set one more independent property to fully specify the state.  For example:

.. code-block:: python

   from sandlersteam.state import State
   state = State(T=200.0, P=0.1)

   print(f"Specific volume: {state.v}") # it is 2.172 m^3/kg

   # Overwrite a property that is not an input
   state.v = 5.5  # this will wipe out previous T and P inputs; state is now underspecified

   # Need to set one more independent property again
   state.P = 0.1  # now state is fully specified again at P=0.1 MPa and v=5.5 m^3/kg

   print(f"New state temperature: {state.T}") # it will be a little over 918.6 C

You can initialize a state with no inputs, then set any two independent properties later:

.. code-block:: python

   from sandlersteam.state import State
   state = State()  # no inputs yet

   state.T = 150.0  # set temperature
   state.P = 0.2    # set pressure, now state is fully specified

   print(f"Specific volume: {state.v}") # it will be 0.9596 m^3/kg

The attributes you can set using dot notation are:

* T : Temperature (°C)
* P : Pressure (MPa)
* v : Specific volume (m³/kg)
* h : Specific enthalpy (kJ/kg)
* s : Specific entropy (kJ/(kg·K))
* u : Specific internal energy (kJ/kg)
* x : Vapor mass fraction (0 to 1 for two-phase states)
* name : just a label for the state (string)

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
   T  =   800 degC                T  =  1200 degC
   P  =     4 MPa                 P  =   4.5 MPa
   v  =  0.12287 m**3 / kg        v  =  0.15098 m**3 / kg
   s  =  7.8502 kJ / (kg * K)     s  =  8.5825 kJ / (kg * K)
   h  =  4141.5 kJ / kg           h  =  5136.9 kJ / kg
   u  =   3650 kJ / kg            u  =  4457.5 kJ / kg
   Pv =  491.48 kJ / kg           Pv =  679.41 kJ / kg

   Property changes:
   ΔT  =   400 degC
   ΔP  =   0.5 MPa
   Δv  =  0.02811 m**3 / kg
   Δs  =  0.7323 kJ / (kg * K)
   Δh  =  995.4 kJ / kg
   Δu  =  807.5 kJ / kg
   ΔPv =  187.93 kJ / kg

From Python
~~~~~~~~~~~

State-change calculations can be performed by creating two State objects representing the initial and final states, then using the built-in :meth:`delta` method to compute property differences:

.. code-block:: python

   from sandlersteam.state import State
   
   # State 1
   state1 = State(T=350, P=7.5)
   
   # State 2
   state2 = State(T=400, P=15.5)
   
   # Calculate property differences using built-in methods
   deltas = state1.delta(state2)
   
   print(f"Δh = {deltas['h']:7g}")
   print(f"Δs = {deltas['s']:7g}")
   print(f"Δu = {deltas['u']:7g}")

Next Steps
----------

* Learn more about the :doc:`cli` for advanced command-line usage
* Explore :doc:`examples` for more complex calculations
* Check the :doc:`api/API` for complete API documentation
