.. _cli:

Command-Line Interface
======================

sandlersteam provides a command-line interface (CLI) for quick thermodynamic calculations without writing Python code.

Overview
--------

The CLI has two main subcommands:

* ``state``: Calculate properties at a single thermodynamic state
* ``delta``: Calculate property changes between two states
* ``avail``: List values of T and P for which explicit steam data is in the tables

Global Options
--------------

.. code-block:: bash

   sandlersteam --help

Shows general help and lists available subcommands. (Note that ``-h`` is reserved for enthalpy inputs.)  

.. code-block:: bash

   sandlersteam --version

Displays the installed version of sandlersteam.

``state`` Command
-----------------

Calculate thermodynamic properties for a pure substance at a specified state.

Syntax
~~~~~~

.. code-block:: bash

   sandlersteam state [OPTIONS]

Required Arguments
~~~~~~~~~~~~~~~~~~

You must provide two of the following state variables:

.. option:: -T, --temperature TEMP

   Temperature (C)

.. option:: -P, --pressure PRESS

   Pressure (MPa)

.. option:: -u, --internal-energy U

   Internal energy (kJ/kg)

.. option:: -v, --specific-volume V

   Specific volume (m³/kg)

.. option:: -h, --enthalpy H

   Enthalpy (kJ/kg)

.. option:: -s, --entropy S

   Entropy (kJ/kg-K)

These are optional arguments:

.. option:: -pu , --pressure-unit PRESSURE_UNIT

   Pressure unit for input and output (default: MPa). Options: MPa, bar, atm

.. option:: -vu, --volume-unit VOLUME_UNIT

   Volume unit for output (default: m3). Options: m3, L, cm3

.. option:: -eu, --energy-unit ENERGY_UNIT

   Energy unit for output (default: kJ). Options: J, kJ, MJ

.. option:: -mu, --mass-unit MASS_UNIT

   Mass unit for output (default: kg). Options: kg, g, lbm, mol

Examples
~~~~~~~~

**Basic calculation:**  Methane at 400 K and 0.5 MPa using the Peng-Robinson EOS:

.. code-block:: bash

   sandlersteam state -T 400 -P 0.5

Output::

   T  =   400 C
   P  =   0.5 MPa
   v  =  0.6173 m3/kg
   s  =  7.7938 kJ/kg-K
   h  =  3271.9 kJ/kg
   u  =  2963.2 kJ/kg
   Pv =  308.65 kJ/kg

**Saturated steam:**  At 0.5 MPa:

.. code-block:: bash

   sandlersteam state -P 0.5 -x 1.0

Output::

   T   =  151.86 C
   P   =   0.5 MPa
   v   =  0.3749 m3/kg
   s   =  6.8213 kJ/kg-K
   h   =  2748.7 kJ/kg
   u   =  2561.2 kJ/kg
   Pv  =  187.45 kJ/kg
   x   = 1.0 kg vapor/kg total
   vL  =  0.001093 m3/kg
   sL  =  1.8607 kJ/kg-K
   hL  =  640.23 kJ/kg
   uL  =  639.68 kJ/kg
   PvL =  0.5465 kJ/kg
   vV  =  0.3749 m3/kg
   sV  =  6.8213 kJ/kg-K
   hV  =  2748.7 kJ/kg
   uV  =  2561.2 kJ/kg
   PvV =  187.45 kJ/kg

``delta`` Command
-----------------

Calculate changes in thermodynamic properties between two states.

Syntax
~~~~~~

.. code-block:: bash

   sandlersteam delta [OPTIONS]

Options
~~~~~~~

You are required to provide state variables for both State 1 and State 2. Use the same options as for the ``state`` command, but append ``1`` or ``2`` to indicate the state.  In addition, 

.. option:: --show-states

   Display full state information for both states in addition to property differences

Examples
~~~~~~~~

**Basic state change calculations:**

.. code-block:: bash

   sandlersteam delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 

Output::

   State-change calculations for water/steam:
   ΔT  =    50 C
   ΔP  =     8 MPa
   Δv  = -0.0175864 m3/kg
   Δs  = -0.33006 kJ/kg-K
   Δh  = -40.67 kJ/kg
   Δu  = -28.99 kJ/kg
   ΔPv = -11.8292 kJ/kg

**With state details:**

.. code-block:: bash

   sandlersteam delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 --show-states

Output::

   State-change calculations for water/steam:

   State 1:                       State 2:
   T  =   350 C                   T  =   400 C
   P  =   7.5 MPa                 P  =  15.5 MPa
   v  =  0.032595 m3/kg           v  =  0.0150086 m3/kg
   s  =  6.1792 kJ/kg-K           s  =  5.84914 kJ/kg-K
   h  =  3001.65 kJ/kg            h  =  2960.98 kJ/kg
   u  =  2758.55 kJ/kg            u  =  2729.56 kJ/kg
   Pv =  244.463 kJ/kg            Pv =  232.633 kJ/kg

   Property changes:
   ΔT  =    50 C
   ΔP  =     8 MPa
   Δv  = -0.0175864 m3/kg
   Δs  = -0.33006 kJ/kg-K
   Δh  = -40.67 kJ/kg
   Δu  = -28.99 kJ/kg
   ΔPv = -11.8292 kJ/kg

Units Reference
---------------

The CLI uses units native to the steam tables in *Chemical Engineering Thermodynamics* by Sandler et al.  Default units are as follows:

=================== ==================
Property            Unit
=================== ==================
Temperature (T)     Celsius (C)
Pressure (P)        Megapascal (MPa)
Specific volume (v) m³/kg
Enthalpy (H)        kJ/kg
Entropy (S)         kJ/(kg-K)
Internal energy (U) kJ/kg
=================== ==================

