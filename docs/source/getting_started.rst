Getting Started
===============

Installation
############

Install via mamba (recommended)
-------------------------------

To install the latest release of imdclient from conda-forge::

    mamba install -c conda-forge imdclient

Install via pip
---------------

To install the latest release of imdclient from PyPI::

    pip install imdclient

Building from Source
####################

We highly recommend using virtual environments to source-build IMDClient. If possible, we recommend that you use `mamba <https://mamba.readthedocs.io/en/latest/>`_ as your package manager through `miniforge <https://github.com/conda-forge/miniforge>`_.

Source build with mamba
-----------------------

Ensure that you have `mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_ installed.

Create and activate a new environment::

    mamba create --name imdclient
    mamba activate imdclient

(Recommended) Install all dependencies using the provided environment YAML files for a clean and reproducible environment::

    mamba env update --name imdclient --file devtools/conda-envs/test_env.yaml
    mamba env update --name imdclient --file docs/requirements.yaml

Build this package from source::

    pip install -e <path/to/repo>

Source build with pip
---------------------

(Optional) Create and activate a virtual environment::

    python -m venv venv
    source venv/bin/activate

Install imdclient from source::

    pip install <path/to/repo>

Development environment installation
------------------------------------
For development or documentation builds, use the following commands after activating your environment:

To install development and documentation dependencies::

    pip install -e <path/to/repo>[doc,test]

To update the development dependencies using mamba::

    mamba env update --name imdclient --file devtools/conda-envs/test_env.yaml

Or to update the documentation building dependencies::

    mamba env update --name imdclient --file docs/requirements.yaml