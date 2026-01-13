.. _contributing:

Contributing to sandlersteam
==============================

Thank you for your interest in contributing to sandlersteam! This document provides guidelines for contributing to the project.

Getting Started
---------------

1. **Fork the repository** on GitHub: https://github.com/cameronabrams/sandlersteam

2. **Clone your fork locally:**

   .. code-block:: bash

      git clone https://github.com/YOUR-USERNAME/sandlersteam.git
      cd sandlersteam

3. **Create a branch** for your changes:

   .. code-block:: bash

      git checkout -b feature/your-feature-name

4. **Install in development mode:**

   .. code-block:: bash

      pip install -e .
      pip install -e ".[dev]"  # Install development dependencies

Development Environment
-----------------------

Required Tools
~~~~~~~~~~~~~~

* Python 3.7 or later
* Git
* A text editor or IDE (VS Code, PyCharm, etc.)

Recommended Setup
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   
   # Install package in editable mode with dev dependencies
   pip install -e ".[dev]"
   
   # Install pre-commit hooks (if configured)
   pre-commit install

Dependencies
~~~~~~~~~~~~

Runtime dependencies:

* numpy
* scipy
* sandlermisc (for for general formatting and utilities)

Development dependencies (optional):

* pytest (testing)
* black (code formatting)
* mypy (type checking)
* sphinx (documentation)
* sphinx_rtd_theme (documentation theme)

Types of Contributions
----------------------

Bug Reports
~~~~~~~~~~~

If you find a bug:

1. Check if it's already reported in `GitHub Issues <https://github.com/cameronabrams/sandlersteam/issues>`_
2. If not, create a new issue with:
   
   * Clear, descriptive title
   * Steps to reproduce
   * Expected behavior
   * Actual behavior
   * Python version and operating system
   * sandlersteam version (``sandlersteam --version``)

Example bug report template::

   **Bug Description**
   Brief description of the bug
   
   **To Reproduce**
   ```python
   from sandlersteam.state import State
   s = State(T=800, P=40).lookup()
   # Steps that trigger the bug
   ```
   
   **Expected Behavior**
   What you expected to happen
   
   **Actual Behavior**
   What actually happened (include error messages)
   
   **Environment**
   - OS: [e.g., Windows 10, Ubuntu 20.04]
   - Python version: [e.g., 3.9.5]
   - sandlersteam version: [e.g., 0.5.0]

Feature Requests
~~~~~~~~~~~~~~~~

To suggest a new feature:

1. Open an issue with the "enhancement" label
2. Describe the feature and its use case
3. Explain why it would be valuable
4. Consider implementation challenges

Code Contributions
~~~~~~~~~~~~~~~~~~

1. Start with an issue (bug or feature) to discuss the change
2. Follow the coding standards below
3. Add tests for new functionality
4. Update documentation
5. Submit a pull request

Coding Standards
----------------

Code Style
~~~~~~~~~~

sandlersteam follows PEP 8 style guidelines:

* Use 4 spaces for indentation (no tabs)
* Maximum line length: 100 characters (flexible for clarity)
* Use descriptive variable names
* Add docstrings to all public functions and classes

Type Hints
~~~~~~~~~~

Use type hints for function signatures:

.. code-block:: python

   from typing import Union, List, Optional
   
   def process_data(
       values: List[float],
       scale: float = 1.0,
       units: Optional[str] = None
   ) -> Union[float, List[float]]:
       """Process numerical data with optional scaling."""
       ...

Docstrings
~~~~~~~~~~

Use NumPy-style docstrings:

.. code-block:: python

   def example_function(param1: int, param2: str) -> bool:
       """
       Brief description of function.
       
       Longer description if needed, explaining what the function does,
       any important considerations, etc.
       
       Parameters
       ----------
       param1 : int
           Description of param1
       param2 : str
           Description of param2
           
       Returns
       -------
       bool
           Description of return value
           
       Raises
       ------
       ValueError
           When param1 is negative
           
       Examples
       --------
       >>> example_function(5, "test")
       True
       
       Notes
       -----
       Any additional notes or warnings.
       
       See Also
       --------
       related_function : Related functionality
       """
       ...

Testing
-------

Running Tests
~~~~~~~~~~~~~

.. code-block:: bash

   # Run all tests
   pytest
   
   # Run specific test file
   pytest tests/test_eos.py
   
   # Run with coverage
   pytest --cov=sandlersteam

Writing Tests
~~~~~~~~~~~~~

Place tests in the ``tests/`` directory:

.. code-block:: python

   # tests/test_new_feature.py
   import pytest
   from sandlersteam.state import State
   
    def test_state_resolve_superheated_exact_at_T_and_P(self):
        state = State(P=0.1, T=300)
        state.lookup()
        self.assertAlmostEqual(state.h, 3074.3, places=1)

        state = State(P=1.0, T=500)
        state.lookup()
        self.assertAlmostEqual(state.s, 7.7622, places=3)

        state = State(P=12.5, T=600)
        state.lookup()
        self.assertAlmostEqual(state.u, 3225.4, places=1)

    def test_state_resolve_superheated_interpolated_at_T_and_P(self):
        state = State(P=0.5, T=375)
        state.lookup()
        self.assertAlmostEqual(state.h, 3219.8, places=1)

        state = State(P=0.7, T=450)
        state.lookup()
        self.assertAlmostEqual(state.s, 7.787225, places=3)

Test Coverage
~~~~~~~~~~~~~

Aim for good test coverage of new code:

* Test normal cases
* Test edge cases
* Test error conditions
* Test with different EOS types

Documentation
-------------

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cd docs
   make html
   # Open docs/_build/html/index.html in browser

Documentation Standards
~~~~~~~~~~~~~~~~~~~~~~~

* Use clear, concise language
* Include code examples
* Cross-reference related functions/classes
* Update changelog for user-facing changes

Adding API Documentation
~~~~~~~~~~~~~~~~~~~~~~~~

API documentation is auto-generated from docstrings. To add new API docs:

1. Write clear docstrings in the source code
2. Add the module/class/function to the appropriate ``docs/api/*.rst`` file
3. Use ``.. autofunction::`` or ``.. autoclass::`` directives

Example:

.. code-block:: rst

   New Function
   ~~~~~~~~~~~~
   
   .. autofunction:: sandlersteam.utils.new_function

Pull Request Process
--------------------

1. **Update your branch:**

   .. code-block:: bash

      git fetch origin
      git rebase origin/main

2. **Run tests and linting:**

   .. code-block:: bash

      pytest
      black sandlersteam tests  # Format code
      mypy sandlersteam  # Type checking (if configured)

3. **Commit your changes:**

   .. code-block:: bash

      git add .
      git commit -m "Brief description of changes"

   Write clear commit messages:
   
   * Use present tense ("Add feature" not "Added feature")
   * First line: brief summary (50 chars or less)
   * Leave blank line, then detailed explanation if needed

4. **Push to your fork:**

   .. code-block:: bash

      git push origin feature/your-feature-name

5. **Create pull request on GitHub:**

   * Clear title describing the change
   * Reference related issues (``Fixes #123``)
   * Describe what changed and why
   * Include test results if applicable
   * Note any breaking changes

Pull Request Template::

   ## Description
   Brief description of changes
   
   ## Related Issues
   Fixes #(issue number)
   
   ## Changes Made
   - Change 1
   - Change 2
   
   ## Testing
   - [ ] Tests pass locally
   - [ ] Added new tests for this feature
   - [ ] Documentation updated
   
   ## Checklist
   - [ ] Code follows project style
   - [ ] Self-reviewed the code
   - [ ] Commented complex code
   - [ ] Updated documentation
   - [ ] No new warnings

Review Process
--------------

What to Expect
~~~~~~~~~~~~~~

* Maintainers will review your PR
* They may request changes or ask questions
* Discussion is encouraged - explain your approach
* Be patient - reviews take time
* Be open to feedback

Responding to Feedback
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Make requested changes
   git add .
   git commit -m "Address review comments"
   git push origin feature/your-feature-name

Code of Conduct
---------------

* Be respectful and inclusive
* Focus on constructive criticism
* Accept feedback gracefully
* Help others learn and improve
* Maintain professionalism

Unacceptable behavior includes:

* Harassment or discrimination
* Trolling or insulting comments
* Personal or political attacks
* Publishing others' private information

Release Process
---------------

For maintainers:

1. Update version in ``pyproject.toml``
2. Update ``docs/changelog.rst``
3. Create git tag: ``git tag v0.X.X``
4. Push tag: ``git push origin v0.X.X``
5. Build and upload to PyPI:

   .. code-block:: bash

      python -m build
      twine upload dist/*

Recognition
-----------

Contributors are recognized in:

* GitHub contributors page
* Release notes for significant contributions
* Special acknowledgments for major features

Questions?
----------

* Open an issue for questions about contributing
* Email: cfa22@drexel.edu
* Check existing issues and pull requests

Thank you for contributing to sandlersteam! Your efforts help make this tool better for everyone in the chemical engineering education community.
