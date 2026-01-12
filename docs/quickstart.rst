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

.. code-block:: bash

   sandlersteam state -TC 400 -P 0.5

Output::

   THERMODYNAMIC STATE OF UNSATURATED STEAM/WATER:
   TC =  400.0 C =  673.1 K
   P  =  0.50 MPa =  5.00 bar
   u  =  2963.2 kJ/kg =  53382.9 J/mol
   v  =  0.6173 m3/kg =  0.0111208 m3/mol
   s  =  7.7938 kJ/kg-K =  140.407 J/mol-K
   h  =  3271.9 kJ/kg =  58944.2 J/mol

From Python
~~~~~~~~~~~

.. code-block:: python

   from sandlersteam.state import State
   state1 = State(TC=200.0, P=0.1)
   print(f"Specific volume: {state1.v} m³/kg")   

You need not set the state variables (T, P) during initialization; you can also just directly assign them later:

.. code-block:: python

   from sandlersteam import PengRobinsonEOS
   # Create EOS object without state and without compound
   pr = PengRobinsonEOS()
   
   # Set compound later; .set_compound returns self for convenience
   pr_methane = pr.set_compound('methane')
   
   # Set state later
   pr_methane.T = 400
   pr_methane.P = 0.5
   
   # Now access properties as before
   print(f"Molar volume: {', '.join([f'{v: 6g}' for v in pr_methane.v])} m³/mol")

State-Change Calculations
-------------------------

Calculate property changes between two thermodynamic states:

From the Command Line
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sandlersteam delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states

Output::

   State-change calculations for methane using Peng-Robinson Equation of State:

   State 1:                       State 2:
   T    =  350.00 K               T    =  400.00 K
   P    =  7.50 MPa               P    =  15.50 MPa
   Z    =  0.92619                Z    =  0.950871
   v    =  0.000359369 m3/mol     v    =  0.000204025 m3/mol
   h    =  929.35 J/mol           h    =  2501.21 J/mol
   s    = -32.095 J/mol-K         s    = -33.5449 J/mol-K
   hdep = -989.935 J/mol          hdep = -1412.32 J/mol
   sdep = -2.12621 J/mol-K        sdep = -2.86197 J/mol-K

   Property changes:
   Δh =  1571.86 J/mol
   Δs = -1.44983 J/mol-K
   Δu =  1104.74 J/mol

From Python
~~~~~~~~~~~

State-change calculations can be performed by creating two EOS objects representing the initial and final states, then using the built-in methods to compute property differences:

.. code-block:: python

   from sandlersteam import PengRobinsonEOS
   
   # State 1
   pr_methane1 = PengRobinsonEOS(T=350, P=7.5).set_compound('methane')
   
   # State 2
   pr_methane2 = PengRobinsonEOS(T=400, P=15.5).set_compound('methane')
   
   # Calculate property differences using built-in methods
   delta_h = pr_methane1.delta_h(pr_methane2)
   delta_s = pr_methane1.delta_s(pr_methane2)
   delta_u = pr_methane1.delta_u(pr_methane2)
   
   print(f"Δh = {', '.join([f'{dh: 7g}' for dh in delta_h])} J/mol")
   print(f"Δs = {', '.join([f'{ds: 7g}' for ds in delta_s])} J/mol-K")
   print(f"Δu = {', '.join([f'{du: 7g}' for du in delta_u])} J/mol")

Available Equations of State
-----------------------------

Ideal gas
~~~~~~~~~~~~~

.. code-block:: python

   from sandlersteam import IdealGasEOS
   
   ig = IdealGasEOS()

van der Waals
~~~~~~~~~~~~~

.. code-block:: python

   from sandlersteam import VanDerWaalsEOS
   
   vdw_methane = VanDerWaalsEOS().set_compound('methane')

Peng-Robinson
~~~~~~~~~~~~~

.. code-block:: python

   from sandlersteam import PengRobinsonEOS
   
   pr_benzene = PengRobinsonEOS().set_compound('benzene')

Soave-Redlich-Kwong
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from sandlersteam import SoaveRedlichKwongEOS
   
   srk_ethanol = SoaveRedlichKwongEOS().set_compound('ethanol')

Next Steps
----------

* Learn more about the :doc:`cli` for advanced command-line usage
* Explore :doc:`examples` for more complex calculations
* Read about the :doc:`theory` behind cubic equations of state
* Check the :doc:`api/API` for complete API documentation
