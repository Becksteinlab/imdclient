# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--
The rules for this file:
  * entries are sorted newest-first.
  * summarize sets of changes - don't reproduce every git log comment here.
  * don't ever delete anything.
  * keep the format consistent:
    * do not use tabs but use spaces for formatting
    * 79 char width
    * YYYY-MM-DD date format (following ISO 8601)
  * accompany each entry with github issue/PR number (Issue #xyz)
-->

## [v0.2.0] - 2025-07-18

**Breaking change**: The removal of the IMDReader and the associated
streaming analysis functionality will *break code* (see issue
#53). The IMDReader is being integrated into MDAnalysis and should be
available in MDAnalysis release 2.10.0.

### Authors
<!-- GitHub usernames of contributors to this release -->
@amruthesht @ljwoods2 @hmacdope @jaclark5 @orbeckst

### Changed
<!-- Changes in existing functionality -->
* IMDReader removed from imdclient by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/54 (issue #53)
* added support for Python 3.13 (issue #90)

### Added
<!-- New added features -->
* Updated installation instructions by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/74
* MD engine links updated in IMDClient documentation by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/73
* Testing: Parse input files for DT when not in traj by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/71

### Fixed
<!-- Bug fixes -->
* updated deployment workflow to use latest pypa/gh-action-pypi-publish@v1.12.4 action (PR #78)


### Deprecated
<!-- Soon-to-be removed features -->

### Removed
<!-- Removed features -->



## [v0.1.5-alpha] - 2025-06-24

### Authors
<!-- GitHub usernames of contributors to this release -->
@orbeckst @amruthesht @hcho38 @ljwoods2

### Added
<!-- New added features -->
* NAMD Testing with private image + Docs by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/46
* Lammps usage docs by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/49
* Instructions for various namd builds by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/51
* chore: added client.py files and cleaned uo output data from examples by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/52
* add codecov token to workflow by @orbeckst in https://github.com/Becksteinlab/imdclient/pull/60
* add sphinx.configuration for RTD by @orbeckst in https://github.com/Becksteinlab/imdclient/pull/62
* Updated CHANGELOG based on the release notes of IMDClient by @hcho38 in https://github.com/Becksteinlab/imdclient/pull/67
* `IMD_TIME` packet definitions modified by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/68

### Changed
<!-- Changes in existing functionality -->
* Run all simulation enegine tests regardless of failed tests by @amruthesht in https://github.com/Becksteinlab/imdclient/pull/69



## [v0.1.4] - 2024-12-13

### Authors
<!-- GitHub usernames of contributors to this release -->
@ljwoods2 @hcho38

### Added
<!-- New added features -->
* Simulation engine GPU + MPI testing by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/21
* Add continue_after_disconnect option (wait packet implementation) by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/44
* Reader test by @hcho38 in https://github.com/Becksteinlab/imdclient/pull/16


## [v0.1.3] - 2024-11-29

### Authors
<!-- GitHub usernames of contributors to this release -->
@ljwoods2 @orbeckst

### Added
<!-- New added features -->
* Docker by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/30

### Fixed
<!-- Bug fixes -->
* Better error messages, Context manager by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/34

### Changed
<!-- Changes in existing functionality -->
* Renamed IMDREADER to IMD by @orbeckst in https://github.com/Becksteinlab/imdclient/pull/38


## [v0.1.2] - 2024-10-28

### Authors
<!-- GitHub usernames of contributors to this release -->
@ljwoods2

### Added
<!-- New added features -->
* Wait flag writeup by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/28

### Fixed
<!-- Bug fixes -->
* URI, Docstring fixes by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/29


## [v0.1.1] - 2024-10-14

### Authors
<!-- GitHub usernames of contributors to this release -->
@ljwoods2

### Changed
<!-- Changes in existing functionality -->
* Changed license from MIT to GPLv3 for compatibility with MDAnalysis by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/25


## [v0.1.0] - 2024-10-13

### Authors
<!-- GitHub usernames of contributors to this release -->
@ljwoods2 @hcho38

### Added
<!-- New added features -->
* New order by @hcho38 in https://github.com/Becksteinlab/imdclient/pull/14
* Manual integration test tool by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/18
* Installation instructions by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/19
* New protocol writeup by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/17
* Stackable, streamable analysis by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/20
* Prepare CI for 0.1.0 release by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/22

### Fixed
<!-- Bug fixes -->
* StreamReaderBase, GROMACS fixes by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/15

### Changed
<!-- Changes in existing functionality -->
* LAMMPS updates by @ljwoods2 in https://github.com/Becksteinlab/imdclient/pull/13
