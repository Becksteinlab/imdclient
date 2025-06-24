IMDClient
==============================
[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  [![Documentation Status][badge_docs]][url_docs]|
| :----------------- | :------- |
| **Status**         | [![GH Actions Status][badge_actions]][url_actions] [![codecov][badge_codecov]][url_codecov] |
| **Community**      | [![License: MIT][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda]|

[badge_actions]: https://github.com/becksteinlab/imdclient/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/becksteinlab/imdclient/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/becksteinlab/imdclient/latest
[badge_docs]: https://readthedocs.org/projects/imdclient/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-MIT-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_release]: https://img.shields.io/github/release-pre/becksteinlab/imdclient.svg
[url_actions]: https://github.com/becksteinlab/imdclient/actions?query=branch%3Amain+workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/becksteinlab/imdclient/branch/main
[url_docs]: https://imdclient.readthedocs.io/en/latest/?badge=latest
[url_latest_release]: https://github.com/becksteinlab/imdclient/releases
[url_license]: https://opensource.org/license/mit
[url_mda]: https://www.mdanalysis.org

Receiver for [IMDv3 protocol](https://imdclient.readthedocs.io/en/latest/protocol_v3.html) from simulation engines like Gromacs, LAMMPS, and NAMD.

IMDClient is bound by a [Code of Conduct](https://github.com/becksteinlab/imdreader/blob/main/CODE_OF_CONDUCT.md).

### Installation

#### Install via conda (recommended)

To install the latest release of IMDClient from conda-forge:

```
conda install -c conda-forge imdclient
```

#### Install via pip

To install the latest release of IMDClient from PyPI:

```
pip install imdclient
```

---

### Building from Source
To build IMDClient from source, we highly recommend using virtual environments.
If possible, we strongly recommend that you use [Anaconda](https://docs.conda.io/en/latest/) as your package manager.

Below we provide instructions both for `conda` and for `pip`.

#### Source build with conda

1. Ensure you have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.

2. Create and activate a new environment:

```
conda create --name imdclient
conda activate imdclient
```

3. Build and install IMDClient in editable mode:

```
pip install -e .
```

4. (Optional) Update dependencies:

```
conda update --all
```

5. Deactivate the environment when finished:

```
conda deactivate
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

Copyright (c) 2024, Lawson


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.

**If you use IMDClient in your research, please cite [IMDClient](https://github.com/Becksteinlab/imdclient) in your publications.**
