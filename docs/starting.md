---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.8.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Getting Started

## Installation

The preferred method to install  `spextra` is using `pip`:

```bash
pip install spextra
```

To install the development version, it is possible to clone the repository and run a local install:

```bash
git clone https://github.com/AstarVienna/speXtra.git
cd speXtra
pip install -e .
```

Keep in mind that the current dev version installed from GitHub may not be stable!
We recommend this only for advanced users, if you know _why_ you need to install the dev version.

### Dependencies

The following dependencies are necessary to run `speXtra`

- [numpy](http://www.numpy.org/)
- [astropy](http://www.astropy.org)
- [scipy](http://www.scipy.org/)
- [synphot](http://synphot.readthedocs.io)
- [PyYAML](https://pyyaml.org/)
- [pooch](https://www.fatiando.org/pooch/)
- [tqdm](https://tqdm.github.io/)
- [astar-utils](https://pypi.org/project/astar-utils/)

Additionally, you may need the following libraries for specific purposes.

- [matplotlib](http://www.matplotlib.org/) for plotting
- [specutils](specutils.readthedocs.io/) optional for loading external spectra in other formats (not required for basic function)

## Basic functionality

The core of `speXtra` is the `Spextrum` which is a wrapper of the `synphot.SourceSpectrum`
with added functionalities.

For example to load a `S0` galaxy templates from the  `kc96` just type

```{code-cell} ipython3
from spextra import Spextrum
sp = Spextrum("kc96/s0")
sp.plot()
```

The last statement will create a plot of the spectrum for a quick visualization

All operations available in {class}`synphot.SourceSpectrum` are possible:

```{code-cell} ipython3
sp1 = Spextrum("kc96/s0")
sp2 = Spextrum("agn/qso")
sp = sp1 + 0.3*sp2
```

The new spectrum will be the sum of the two components

### Scaling to a magnitude

```{code-cell} ipython3
from astropy import units as u
sp1 = Spextrum("kc96/s0")
sp2 = sp1.scale_to_magnitude(amplitude=13 * u.ABmag, filter_curve="g")
```

### Obtaining magnitudes from spectra

```{code-cell} ipython3
mag = sp1.get_magnitude(filter_curve="g")
print(mag)
```

### Redshifting the spectra

It is possible to specify a redshift `z` ir a velocity `vel` (negative velocities are allowed)

```{code-cell} ipython3
sp3 = sp2.redshift(z=1)
vel = -1000 * u.km / u.s
sp2 = sp1.redshift(vel=vel)
```
