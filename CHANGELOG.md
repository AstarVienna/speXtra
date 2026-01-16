# Version 0.43.0
**2026-01-16**

Add support for Python 3.14, minor bugfixes.

## What's Changed
### Bugs fixed
* Don't include "*_unit" keywords for FITS spectra by @teutoburg in https://github.com/AstarVienna/speXtra/pull/97
### Dependency Changes
* Bump fonttools from 4.47.2 to 4.61.0 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/95
* Bump urllib3 from 2.5.0 to 2.6.0 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/96
* Bump urllib3 from 2.6.0 to 2.6.3 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/98
* Allow Python 3.14, sync dependencies by @teutoburg in https://github.com/AstarVienna/speXtra/pull/99
### Other Changes
* Minor cleanup in tests by @teutoburg in https://github.com/AstarVienna/speXtra/pull/100
* Turn `DeprecationWarning`s into `FutureWarning`s, specify version by @teutoburg in https://github.com/AstarVienna/speXtra/pull/101

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.42.0...v0.43.0


# Version 0.42.0
**2025-09-17**

Add support for Python 3.13 and Numpy 2, solve some other dependency issues.

## What's Changed
### Dependency Changes
* Bump requests from 2.32.2 to 2.32.4 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/79
* Bump urllib3 from 2.2.2 to 2.5.0 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/80
* Bump actions/checkout from 4 to 5 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/85
* Bump actions/download-artifact from 4.3.0 to 5.0.0 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/84
* Migrate to Poetry 2, allow numpy 2.x by @teutoburg in https://github.com/AstarVienna/speXtra/pull/87
* Bump actions/setup-python from 5 to 6 by @dependabot[bot] in https://github.com/AstarVienna/speXtra/pull/89
* Add support for Python 3.13 by @teutoburg in https://github.com/AstarVienna/speXtra/pull/92
### Documentation Improvements
* Fix readme and toml by @teutoburg in https://github.com/AstarVienna/speXtra/pull/91
### Other Changes
* Fix CI comment by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/83

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.41.4...v0.42.0


# Version 0.41.4
**2025-05-27**

Hotfix to solve some unit issues.

## What's Changed
### Bugs fixed
* Catch dimensionless quantities in wave and flux by @teutoburg in https://github.com/AstarVienna/speXtra/pull/78
### Dependency Changes
* Bump tornado from 6.4.1 to 6.4.2 by @dependabot in https://github.com/AstarVienna/speXtra/pull/65
* Bump actions/download-artifact from 4.1.7 to 4.1.8 by @dependabot in https://github.com/AstarVienna/speXtra/pull/68
* Bump abatilo/actions-poetry from 3 to 4 by @dependabot in https://github.com/AstarVienna/speXtra/pull/69
* Bump actions/download-artifact from 4.1.8 to 4.1.9 by @dependabot in https://github.com/AstarVienna/speXtra/pull/70
* Bump jinja2 from 3.1.5 to 3.1.6 by @dependabot in https://github.com/AstarVienna/speXtra/pull/71
* Bump actions/download-artifact from 4.1.9 to 4.2.0 by @dependabot in https://github.com/AstarVienna/speXtra/pull/72
* Bump actions/download-artifact from 4.2.0 to 4.2.1 by @dependabot in https://github.com/AstarVienna/speXtra/pull/73
* Bump pytest and pytest-cov to match ScopeSim by @teutoburg in https://github.com/AstarVienna/speXtra/pull/74
* Bump actions/download-artifact from 4.2.1 to 4.3.0 by @dependabot in https://github.com/AstarVienna/speXtra/pull/76
* Bump tornado to 6.5 to solve security alert by @teutoburg in https://github.com/AstarVienna/speXtra/pull/77
### Other Changes
* Treat warnings as errors in testing by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/63
* Ignore our own DeprecationWarnings by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/64
* Bump jinja2 by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/67
* Add dependabot configuration by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/66
* Initialize self.repr to aid debugging by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/75

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.41.3...v0.41.4


# Version 0.41.3
**2024-11-01**

Hotfix to solve some combatibility issues.

## What's Changed
### Dependency Changes
* Bump some dependencies to match ScopeSim by @teutoburg in https://github.com/AstarVienna/speXtra/pull/58
### Other Changes
* Fix waves check by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/57
* Add doctests to workflow by @teutoburg in https://github.com/AstarVienna/speXtra/pull/59

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.41.2...v0.41.3


# Version 0.41.2
**2024-10-26**

## What's Changed
### Dependency Changes
* Update to astropy 6.0.1 by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/54
### Other Changes
* Allow dashes in filter names by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/52
* Fix webtests by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/55
* Run jobs nightly by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/53

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.41.1...v0.41.2


# Version 0.41.1
**2024-09-09**

Hotfix for v0.41.0 to include images in README on PyPI plus some internal cleanup.

## What's Changed
### Documentation Improvements
* Include static docs images in dist to fix readme by @teutoburg in https://github.com/AstarVienna/speXtra/pull/47
### Other Changes
* Remove two obsolete files by @teutoburg in https://github.com/AstarVienna/speXtra/pull/48
* Improve publish workflow by @teutoburg in https://github.com/AstarVienna/speXtra/pull/50

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.41.0...v0.41.1

# Version 0.41.0
**2024-09-09**

Mostly internal maintenance, and some functions are now stricter in requiring quantified input parameters (aka astropy units).

> [!IMPORTANT]
> The minimum required Python version for this package is now **3.10** (see Dependency Changes).

## What's Changed
### API Changes
* Quantities and more by @teutoburg in https://github.com/AstarVienna/speXtra/pull/43
### Dependency Changes
* Bump pillow from 10.2.0 to 10.3.0 by @dependabot in https://github.com/AstarVienna/speXtra/pull/26
* Bump idna from 3.6 to 3.7 by @dependabot in https://github.com/AstarVienna/speXtra/pull/27
* Bump jinja2 from 3.1.3 to 3.1.4 by @dependabot in https://github.com/AstarVienna/speXtra/pull/29
* Bump tqdm from 4.66.1 to 4.66.3 by @dependabot in https://github.com/AstarVienna/speXtra/pull/28
* Bump astropy, requests, idna, pillow by @teutoburg in https://github.com/AstarVienna/speXtra/pull/31
* Bump tornado from 6.4 to 6.4.1 by @dependabot in https://github.com/AstarVienna/speXtra/pull/32
* Bump urllib3 from 2.1.0 to 2.2.2 by @dependabot in https://github.com/AstarVienna/speXtra/pull/33
* Bump certifi from 2023.11.17 to 2024.7.4 by @dependabot in https://github.com/AstarVienna/speXtra/pull/34
* Bump zipp from 3.17.0 to 3.19.1 by @dependabot in https://github.com/AstarVienna/speXtra/pull/35
* Use `more_itertools.always_iterable` instead of DIY function by @teutoburg in https://github.com/AstarVienna/speXtra/pull/38
* Bump actions/download-artifact from 3 to 4.1.7 in /.github/workflows by @dependabot in https://github.com/AstarVienna/speXtra/pull/40
* Drop support for Python 3.9 by @teutoburg in https://github.com/AstarVienna/speXtra/pull/41
### Other Changes
* Use new main branch for DevOps workflows by @teutoburg in https://github.com/AstarVienna/speXtra/pull/36
* Remove remaining master references by @teutoburg in https://github.com/AstarVienna/speXtra/pull/39
* Update workflow versions in publish_pypi by @teutoburg in https://github.com/AstarVienna/speXtra/pull/42
* Call bump workflow in devops to avoid copies by @teutoburg in https://github.com/AstarVienna/speXtra/pull/44

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.40.0...v0.41.0

# Version 0.40.0
**2024-03-19**

This version introduces a radical restructuring of the package's internal structure regarding the object model and MRO of the "database" part that allows creation of the main classes from library references. The public API is mostly unchanged, but there are several deprecations to look out for. The documentation is mostly updated for these changes, but might still include some references to the old structure.

## What's Changed
### API Changes
* Complete Restructuring by @teutoburg in https://github.com/AstarVienna/speXtra/pull/10
### Dependency Changes
* Bump jinja2 from 3.1.2 to 3.1.3 by @dependabot in https://github.com/AstarVienna/speXtra/pull/12
* Drop support for Python 3.8 by @teutoburg in https://github.com/AstarVienna/speXtra/pull/13
* Allow internal updates by @teutoburg in https://github.com/AstarVienna/speXtra/pull/17
### Documentation Improvements
* Add config file for auto-generated release notes by @teutoburg in https://github.com/AstarVienna/speXtra/pull/11
* Multiple improvements to packaging and release by @teutoburg in https://github.com/AstarVienna/speXtra/pull/14
* Fix notebooks so they run again by @teutoburg in https://github.com/AstarVienna/speXtra/pull/16
* Update docs for RTD by @teutoburg in https://github.com/AstarVienna/speXtra/pull/18
* Cleanup and improve documentation on readthedocs by @teutoburg in https://github.com/AstarVienna/speXtra/pull/21
* Add sphinx-copybutton to docs by @teutoburg in https://github.com/AstarVienna/speXtra/pull/22
* Include pdf and epub in docs build by @teutoburg in https://github.com/AstarVienna/speXtra/pull/23
* Fix RTD poetry configuration by @teutoburg in https://github.com/AstarVienna/speXtra/pull/24
### Other Changes
* Use pyproject.toml and DevOps workflow by @hugobuddel in https://github.com/AstarVienna/speXtra/pull/9
* Allow most tests to run offline by @teutoburg in https://github.com/AstarVienna/speXtra/pull/15
* Use split tests from DevOps by @teutoburg in https://github.com/AstarVienna/speXtra/pull/19
* Fix emission line spectrum by @teutoburg in https://github.com/AstarVienna/speXtra/pull/25

## New Contributors
* @teutoburg made their first contribution in https://github.com/AstarVienna/speXtra/pull/11
* @dependabot made their first contribution in https://github.com/AstarVienna/speXtra/pull/12

**Full Changelog**: https://github.com/AstarVienna/speXtra/compare/v0.33...v0.40.0

# Version 0.33.0

## What's Changed
- Remove vegamag as a unit.

# Version 0.23.0

Update versioning scheme.

## What's Changed
- Bug in path resolution in windows solved.

# Version 0.2.0

## What's Changed
- Most functionally is there and tested
- tynt requirement dropped
- Docs updated

# Version 0.1.0

Basic functionality is there.
