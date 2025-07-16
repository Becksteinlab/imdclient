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

IMDClient requires Python 3.10 or higher.

#### Install via mamba (recommended)
To install the latest release of IMDClient from conda-forge:

```
mamba install -c conda-forge imdclient
```

#### Install via pip
To install the latest release of IMDClient from PyPI:

```
pip install imdclient
```

---

### Building from Source
To build IMDClient from source, we highly recommend using virtual environments.
If possible, we recommend that you use [mamba](https://mamba.readthedocs.io/en/latest/) as your package manager through [miniforge](https://github.com/conda-forge/miniforge). (You can substitute `conda` for `mamba` in the commands below if you prefer.)

Below we provide instructions both for `mamba` and for `pip`.

#### Source build with mamba
1. Ensure you have [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) installed.

2. Create and activate a new environment:

```
mamba create --name imdclient
mamba activate imdclient
```

3. (Recommended) Install all dependencies using the provided environment YAML files for a clean and reproducible environment:

```
mamba env update --name imdclient --file devtools/conda-envs/test_env.yaml
mamba env update --name imdclient --file docs/requirements.yaml
```

4. Build and install IMDClient in editable mode:

```
pip install -e .
```

5. (Optional) Update dependencies:

```
mamba update --all
```

6. Deactivate the environment when finished:

```
mamba deactivate
```

#### Source build with pip

1. (Optional) Create and activate a virtual environment:

```
python -m venv venv
source venv/bin/activate
```

2. Install IMDClient from source:

```
pip install .
```

3. (Optional) For development (tests and docs):

```
pip install ".[test,doc]"
```

### Copyright

The IMDClient source code is hosted at https://github.com/becksteinlab/imdclient
and is available under the MIT license (see the file [LICENSE](https://github.com/becksteinlab/imdclient/blob/main/LICENSE)).

Copyright (c) 2024-2025, imdclient [AUTHORS](https://github.com/Becksteinlab/imdclient/blob/main/AUTHORS.md)


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.

**If you use IMDClient in your research, please cite [IMDClient](https://github.com/Becksteinlab/imdclient) in your publications.**
