IMDClient
==============================
[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  [![Documentation Status][badge_docs]][url_docs]|
| :----------------- | :------- |
| **Status**         | [![GH Actions Status][badge_actions]][url_actions] [![codecov][badge_codecov]][url_codecov] |
| **Community**      | [![License: MIT][badge_license]][url_license]

[badge_actions]: https://github.com/becksteinlab/imdclient/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/becksteinlab/imdclient/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/becksteinlab/imdclient/latest
[badge_docs]: https://readthedocs.org/projects/imdclient/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-MIT-blue.svg
[badge_release]: https://img.shields.io/github/release-pre/becksteinlab/imdclient.svg
[url_actions]: https://github.com/becksteinlab/imdclient/actions?query=branch%3Amain+workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/becksteinlab/imdclient/branch/main
[url_docs]: https://imdclient.readthedocs.io/en/latest/?badge=latest
[url_latest_release]: https://github.com/becksteinlab/imdclient/releases
[url_license]: https://opensource.org/license/mit

Receiver for [IMDv3 protocol](https://imdclient.readthedocs.io/en/latest/protocol_v3.html) from simulation engines like Gromacs, LAMMPS, and NAMD.

IMDClient is bound by a [Code of Conduct](https://github.com/becksteinlab/imdclient/blob/main/CODE_OF_CONDUCT.md).

### Installation

IMDClient is available via PyPi and can be installed with pip:
```bash
pip install imdclient
```

To build IMDClient from source,
we highly recommend using virtual environments.
If possible, we strongly recommend that you use
[Anaconda](https://docs.conda.io/en/latest/) as your package manager.
Below we provide instructions both for `conda` and
for `pip`.

#### With conda

Ensure that you have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.

Create a virtual environment and activate it:

```
conda create --name imdclient
conda activate imdclient
```

<!-- Install the development and documentation dependencies:

```
conda env update --name imdclient --file devtools/conda-envs/test_env.yaml
conda env update --name imdclient --file docs/requirements.yaml
``` -->

Build this package from source:

```
pip install -e .
```

If you want to update your dependencies (which can be risky!), run:

```
conda update --all
```

And when you are finished, you can exit the virtual environment with:

```
conda deactivate
```

#### With pip

To build the package from source, run:

```
pip install .
```

If you want to create a development environment, install
the dependencies required for tests and docs with:

```
pip install ".[test,doc]"
```

### Copyright

The IMDClient source code is hosted at https://github.com/becksteinlab/imdclient
and is available under the MIT license (see the file [LICENSE](https://github.com/becksteinlab/imdclient/blob/main/LICENSE)).

Copyright (c) 2024, Lawson


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
