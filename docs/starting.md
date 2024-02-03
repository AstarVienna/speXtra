# Getting Started

The core of `speXtra` is the {class}`Spextrum` which is a wrapper of the {class}`synphot.SourceSpectrum`
with added functionalities.

For example to load a `S0` galaxy templates from the  {ref}`kc96` (`kc96`) just type

```python
from spextra import Spextrum
sp = Spextrum("kc96/s0")
sp.plot()
```

The last statement will create a plot of the spectrum for a quick visualization

All operations available in   {class}`synphot.SourceSpectrum` are possible:

```python
sp1 = Spextrum("kc96/s0")
sp2 = Spextrum("agn/qso")
sp = sp1 + 0.3*sp2
```

The new spectrum will be the sum of the two components

## Scaling to a magnitude

```python
sp1 = Spextrum("kc96/s0")
sp2 = sp1.scale_to_magnitude(amplitude=13 * u.ABmag, filter_curve="g")
```

## Obtaining magnitudes from spectra

```python
mag = sp1.get_magnitude(filter_curve="g")
print(mag)
```

## Redshifting the spectra

It is possible to specify a redshift `z` ir a velocity `vel` (negative velocities are allowed)

```python
sp3 = sp2.redshift(z=1)

import astropy.units as u

vel = -1000 * u.km / u.s
sp2 = sp1.redshift(vel=vel)
```
