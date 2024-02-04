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

# speXtra workflow

## Getting Started

Create a spectra using the built-in `'s0'` galaxy source from the Kinney-Calzetti library
`'kc96'`:

```{code-cell} ipython3
  from spextra import Spextrum
  sp = Spextrum("kc96/s0")
```

and `spextra` will download the spectra and create an object

## Creating a spextrum from a file
