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

# Filters

Most of the filters available to `speXtra` are downloaded from
the `SVO database <http://svo2.cab.inta-csic.es/theory/fps/>`. However, many other
system are not included in that database so a few (at the moment) are included in the database.

To list the filter systems included in the database

```{code-cell} ipython3
    from spextra import spextra_database
    spextra_database["filter_systems"]
```

The ::class::`Pasband` holds the information of an astronomical filter and can be called simply using:

```{code-cell} ipython3
from spextra import Passband
filt = Passband("elt/micado/Y")
filt.plot()
```

Many standard filters have been defined using shortcuts for quicker access.
Here is a list.

```{code-cell} ipython3
from spextra import DEFAULT_DATA
DEFAULT_DATA.filters
```

These filters can be used simply using the shortcut

```{code-cell} ipython3
from spextra import Passband

filt = Passband("g")  # for sdss g-band filter
filt.plot()
```

All methods in ::class::`Spextrum` that require to provide a filter name, can be used in either way,
using the full name or just a shortcut.

For a full list of the filters available to `spextra` please visit the [SVO](http://svo2.cab.inta-csic.es/theory/fps/)

Additionally below you can find the filter systems for prospective ELT instruments and others not
available at the SVO

## MICADO

```{literalinclude} ../database/filter_systems/elt/micado/index.yml
:language: yaml
```

## METIS

```{literalinclude} ../database/filter_systems/elt/metis/index.yml
:language: yaml
```

## ETC

Also the standard filters used by the ESO Exposure Time Calculator

```{literalinclude} ../database/filter_systems/etc/index.yml
:language: yaml
```
