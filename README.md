<img src="https://raw.githubusercontent.com/AstarVienna/spextra/main/docs/_static/images/speXtra_logo.png" width="400pt">

# SpeXtra

[![Tests](https://github.com/AstarVienna/speXtra/actions/workflows/tests.yml/badge.svg)](https://github.com/AstarVienna/speXtra/actions/workflows/tests.yml)
[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
![dev version](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FAstarVienna%2FspeXtra%2Fmaster%2Fpyproject.toml&query=%24.tool.poetry.version&label=dev%20version&color=teal)

[![Documentation Status](https://readthedocs.org/projects/spextra/badge/?version=latest)](https://speXtra.readthedocs.io/en/latest)
[![codecov](https://codecov.io/gh/AstarVienna/speXtra/graph/badge.svg)](https://codecov.io/gh/AstarVienna/speXtra)
[![PyPI - Version](https://img.shields.io/pypi/v/speXtra)](https://pypi.org/project/speXtra/)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/speXtra)

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

A python tool to manage and manipulate astronomical spectra


## Description

``speXtra`` is a python tool to download, load, display and manipulate spectra of astronomical sources.
It has developed to provide spectral sources to [ScopeSim](https://scopesim.readthedocs.io/en/latest/) but it may be helpful for other purposes too.

speXtra stands in the shoulder of giants: [synphot](https://synphot.readthedocs.io/en/latest/) and [astropy](https://www.astropy.org/).

To install ``spextra`` simply run:

```bash
pip install spextra
```

## Functionalities

``speXtra`` is able to:

- Download spectra from a database and return it in format compatible with ``synphot`` format.

```python
from spextra import Spextrum
sp = Spextrum("kc96/s0")
```

and it will load the S0 galaxy template of the Kinney-Calzetti spectral library. To quickly
plot the resulting spectra, simply type

```python
sp.plot()
```
<img src="https://raw.githubusercontent.com/AstarVienna/spextra/main/docs/_static/images/kc96_S0.png" width="400pt">
