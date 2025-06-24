Getting Started
===============

Installation
############

Install via conda (recommended)
-------------------------------

To install the latest release of imdclient from conda-forge::

    conda install -c conda-forge imdclient

Install via pip
---------------

To install the latest release of imdclient from PyPI::

    pip install imdclient

Building from Source
####################

We highly recommend using virtual environments to source-build IMDClient. If possible, we strongly recommend that you use `Anaconda <https://docs.conda.io/en/latest/>`_ as your package manager.

Source build with conda
-----------------------

Ensure you have `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ installed.

Create and activate a new environment::

    conda create --name imdclient
    conda activate imdclient

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

To update the development dependencies using conda::

    conda env update --name imdclient --file devtools/conda-envs/test_env.yaml

Or to update the documentation building dependencies::

    conda env update --name imdclient --file docs/requirements.yaml