Getting Started
===============

Installation
############

To build imdclient from source, we highly recommend using virtual environments.
If possible, we strongly recommend that you use
`Anaconda <https://docs.conda.io/en/latest/>`_ as your package manager.
Below we provide instructions both for installing into a `conda` environment.

With conda
----------

Ensure that you have `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ installed.

Create a virtual environment and activate it::

    conda create --name imdclient
    conda activate imdclient

Build this package from source::

    pip install -e <path/to/repo>

Development environment installation
------------------------------------

After creating and activating a conda environment as described, install 
the package with documentation and testing dependencies::

    pip install -e <path/to/repo>[doc, test]

Then, to install the development dependencies::

    conda env update --name imdclient --file devtools/conda-envs/test_env.yaml

Or the documentation building dependencies::

    conda env update --name imdclient --file docs/requirements.yaml