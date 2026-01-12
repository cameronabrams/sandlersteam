.. _installation:

Installation
============

Requirements
------------

sandlersteam requires:

* Python 3.7 or later
* NumPy
* SciPy

Install from PyPI
-----------------

The easiest way to install sandlersteam is using pip:

.. code-block:: bash

   pip install sandlersteam

This will automatically install all required dependencies.

Install from Source
-------------------

To install the latest development version from GitHub:

.. code-block:: bash

   git clone https://github.com/cameronabrams/sandlersteam.git
   cd sandlersteam
   pip install -e .

This installs the package in "editable" mode, which is useful for development.

Verify Installation
-------------------

To verify that sandlersteam is installed correctly, run:

.. code-block:: bash

   sandlersteam --help

You should see the help message for the command-line interface.

Alternatively, from Python:

.. code-block:: python

   import sandlersteam
   print(sandlersteam.__version__)

Updating
--------

To update to the latest version:

.. code-block:: bash

   pip install --upgrade sandlersteam

Uninstalling
------------

To remove sandlersteam:

.. code-block:: bash

   pip uninstall sandlersteam
