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

# Spectral Libraries

A number of spectral libraries have been include in `speXtra`.
See {ref}`database-contents` for a summary what is included in the database.

To list the names of the libraries included in the database

```{code-cell} ipython3
from spextra import spextra_database
spextra_database["libraries"]
```

To see which templates are available in each library

```{code-cell} ipython3
from spextra import SpecLibrary
lib = SpecLibrary("kc96")
list(lib)
```

Below you can find a detailed description of each library.

## A Library of Reference Stars

```{literalinclude} ../database/libraries/ref/index.yml
:language: yaml
```

## The Kinney-Calzetti Spectral Atlas of Galaxies

```{literalinclude} ../database/libraries/kc96/index.yml
:language: yaml
```

## Pickles Stellar Library

```{literalinclude} ../database/libraries/pickles/index.yml
:language: yaml
```

## SDSS galaxy composite spectra

```{literalinclude} ../database/libraries/dobos/index.yml
:language: yaml
```

## IRTF spectral library

```{literalinclude} ../database/libraries/irtf/index.yml
:language: yaml
```

## AGN templates

```{literalinclude} ../database/libraries/agn/index.yml
:language: yaml
```

## Emission line nebulae

```{literalinclude} ../database/libraries/nebulae/index.yml
:language: yaml
```

## Galaxy SEDs from the UV to the Mid-IR

```{literalinclude} ../database/libraries/brown/index.yml
:language: yaml
```

## Subset of Kurucz 1993 Models (subset)

```{literalinclude} ../database/libraries/kurucz/index.yml
:language: yaml
```

## Supernova Legacy Survey

```{literalinclude} ../database/libraries/sne/index.yml
:language: yaml
```

## Flux/Telluric standards with X-Shooter

```{literalinclude} ../database/libraries/moehler/index.yml
:language: yaml
```

## High-Resolution Spectra of Habitable Zone Planets (example)

```{literalinclude} ../database/libraries/madden/index.yml
:language: yaml
```

## BOSZ Stellar Atmosphere Grid (subset) - High Resolution

```{literalinclude} ../database/libraries/bosz/hr/index.yml
:language: yaml
```

## BOSZ Stellar Atmosphere Grid (subset) - Medium Resolution

```{literalinclude} ../database/libraries/bosz/mr/index.yml
:language: yaml
```

## BOSZ Stellar Atmosphere Grid (subset) - Low Resolution

```{literalinclude} ../database/libraries/bosz/lr/index.yml
:language: yaml
```

## Paranal Night Sky Spectra

Additionally, the Paranal sky emission spectra is also included

```{literalinclude} ../database/libraries/sky/index.yml
:language: yaml
```
