dazsle-hbb-recipes
==========================

.. inclusion-marker-1-do-not-remove

This is a highly idiomatic package that contains all the recipes for IDs, generator tools, k-factors, and various scalings for the DAZSLE Hbb/Hcc analysis.

.. inclusion-marker-1-5-do-not-remove

This package is currently organized into five subpackages:

1) gentools - for iterating through baconprod particle histories
2) ids - POG Ids and a few veto selections on top of baconprod
3) kfactors - kfactors for V+jets and TTbar on top of BaconProd
4) plots - plots we want to make in various signal regions
5) weights - scale factors and related

.. inclusion-marker-2-do-not-remove

Installation
============

Install uproot-methods like any other Python package:

.. code-block:: bash

    pip install fnal-column-analysis-tools

or similar (use ``sudo``, ``--user``, ``virtualenv``, or pip-in-conda if you wish).

Strict dependencies:
====================

- `Python <http://docs.python-guide.org/en/latest/starting/installation/>`__ (2.7+, 3.6+)

The following are installed automatically when you install uproot with pip:

- `numpy <https://scipy.org/install.html>`__ (1.15+)
- `awkward-array <https://pypi.org/project/awkward>`__ to manipulate data from non-flat TTrees, such as jagged arrays (`part of Scikit-HEP <https://github.com/scikit-hep/awkward-array>`__)
- `uproot-methods <https://pypi.org/project/uproot-methods>`__ to allow expressions of things as lorentz vectors
- `numba <https://numba.pydata.org/>`__ just-in-time compilation of python functions
- ``scipy`` for statistical functions
- ``matplitlib`` as a plotting backend
- ``uproot`` for interacting with ROOT files

.. inclusion-marker-3-do-not-remove

Tutorial
========

This library is installed by people doing collider HEP analysis in the FNAL CMS group (so far).

Reference documentation
=======================

(...)
