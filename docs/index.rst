.. sandlersteam documentation master file

sandlersteam
=============

.. image:: https://img.shields.io/pypi/v/sandlersteam.svg
   :target: https://pypi.org/project/sandlersteam/
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/pyversions/sandlersteam.svg
   :target: https://pypi.org/project/sandlersteam/
   :alt: Python versions

Digitized steam tables from Sandler's 5th edition
-------------------------------------------------

**sandlersteam** implements a Python interface to the steam tables found in Appendix III of 
*Chemical, Biochemical, and Engineering Thermodynamics* (5th edition) by Stan Sandler (Wiley, USA). 

.. warning::
   This package should be used for **educational purposes only**.

Features
--------

* Command-line interface for quick calculations
* Python API for integration into larger workflows
* Reports key thermodynamic properties of water and steam, including:
  
  * Specific volume (v), 
  * specific enthalpy (h),
  * specific entropy (s), and 
  * specific internal energy (u)

Quick Start
-----------

Installation::

   pip install sandlersteam

Basic usage from the command line::

   sandlersteam state -T 800 -P 40

Basic usage from Python:

.. code-block:: python

   from sandlersteam.state import State
   state = State(T=800, P=40).lookup()
   print(f"Specific volume: {state.v} m³/kg")

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   installation
   quickstart
   cli
   examples

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api/API.rst

.. toctree::
   :maxdepth: 1
   :caption: Developer Guide
   
   contributing
   changelog

License
=======

This project is licensed under the MIT License - see the LICENSE file for details.

Citation
========

If you use this package for academic work, please cite:

* Sandler, S. (2017). *Chemical, Biochemical, and Engineering Thermodynamics* (5th ed.). Wiley.

Contact
=======

Cameron F. Abrams - cfa22@drexel.edu

GitHub: https://github.com/cameronabrams/sandlersteam
