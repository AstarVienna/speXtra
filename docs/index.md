# speXtra

## speXtra

`speXtra` is a `python` library for managing and manipulating astronomical spectra.

Under the hood it packages many `synphot` workflows to perform typical
treatment of astronomical spectra, including scaling, binning, smoothing, redshifting, etc.
`speXtra` does not make measurement on the spectra, except obtaining magnitudes.

`speXtra` comes with a built-in `database` to facilitate the access
to commonly used spectral templates, extinction curves and filters.

:::{warning}
July 2022: The downloadable content server was retired and the data migrated to a new server.

SpeXtra v0.25 and above have been redirected to a new server URL.

For older versions, please either upgrade to the latest version (`pip install --upgrade spextra`), or follow these [instructions to update the server URL](https://astarvienna.github.io/server_upgrade_instructions.html) in the config file.
:::

```{toctree}
:maxdepth: 1

install
starting
notebooks/Database.ipynb
notebooks/speXtra_demo.ipynb
spectral_libraries
extinction_curves
filters
Reference API <reference/spextra>
```
