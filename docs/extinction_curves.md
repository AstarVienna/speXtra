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

# Extinction Curves

Few extinction curves have been included in `speXtra`

To list the names of the extinction curves included in the database

```{code-cell} ipython3
    from spextra import spextra_database
    spextra_database["extinction_curves"]
```

## Gordon

```{literalinclude} ../database/extinction_curves/gordon/index.yml
:language: yaml
```

## Cardelli

```{literalinclude} ../database/extinction_curves/cardelli/index.yml
:language: yaml
```

## Calzetti

```{literalinclude} ../database/extinction_curves/calzetti/index.yml
:language: yaml
```
