# -*- coding: utf-8 -*-
"""Utility functions for downloads.

Cache locations are resolved through ``astar_utils.cache_dir``, which knows
about ScopeSim_Data and the shared ``~/.astar`` cache. All pooch-specific code
lives here in spextra, so that astar-utils stays dependency-free.

The pattern throughout is "check trusted local copies first, download only as a
last resort". Trusted local copies live in ScopeSim_Data (if installed) and in
the package's own bundled data (present in a clone / editable install). The
home cache below them is managed by pooch for registry files, and written to
directly for SVO files.

Registry database
-----------------
Via the module-level ``retriever`` using known hashes. Use the
``fetch_data_file`` wrapper: it pre-checks the trusted local locations and
otherwise delegates to pooch, which owns the home cache (with hash + update
detection) and downloads from the remote database if needed.

SVO filter curves
-----------------
Via ``download_svo_filter``. These have no known hash and are treated as
static: once a copy exists anywhere it is trusted. For now the actual download
uses ``pooch.retrieve(known_hash=None)``. This is a candidate to be replaced by
a reusable ecosystem-wide chunked downloader (e.g. niquests-based) later.
"""

import warnings
from pathlib import Path

import pooch

from astar_utils.cache_dir import find_cached_file

from .configuration import config, __data_dir__

__all__ = [
    "retriever",
    "fetch_data_file",
    "download_file",
    "download_svo_filter",
]



retriever = pooch.create(
    path=config.cache_dir,
    base_url=config.database_url,
    retry_if_failed=config.retry,
)
retriever.load_registry(config.registry_file)


def fetch_data_file(filename: str, **kwargs) -> Path:
    """Return a registry data file, preferring trusted local copies.

    Checks ScopeSim_Data and the package's bundled data first (committed,
    trusted locations, returned as-is without re-downloading). Failing that,
    delegates to pooch, which manages the home download cache -- including hash
    verification and update detection -- and downloads from the remote database
    if the file is missing or stale.

    The home cache is deliberately not pre-checked here: it is pooch's
    ``path``, so pooch should own it (hence ``include_home_cache=False``).
    """
    try:
        cached_path = find_cached_file(
            Path(filename),
            config.package_name,
            extra_dirs=[__data_dir__],
            include_home_cache=False,
        )
        return cached_path
    except FileNotFoundError:
        pass  # Fallback to download
    return Path(retriever.fetch(filename, **kwargs))


def download_file(remote_url, local_name):
    """For backwards compatibility only."""
    warnings.warn("The download_file function is deprecated and will be "
                  "removed in v1.0. Please use retriever.fetch instead.",
                  FutureWarning, stacklevel=2)
    file = pooch.retrieve(remote_url, known_hash=None,
                          fname=local_name.name, path=local_name.parent)
    return file


def download_svo_filter(filter_name):
    """
    Query the SVO service for the transmittance of a given filter.

    The filter is looked for, in order, in ScopeSim_Data, the package's bundled
    data, and the local download cache. Only if it is found nowhere is it
    downloaded from SVO and stored in the write cache. SVO curves are static,
    so an existing copy is always trusted.

    Parameters
    ----------
    filter_name : str
        Name of the filter as available on the Spanish VO filter service,
        e.g. ``"Paranal/HAWKI.Ks"``.

    Returns
    -------
    wave, trans
        Wavelength and transmission columns of the filter curve.
    """
    from astropy.table import Table

    # Relative path that is correct on every platform, e.g.
    # "svo/Paranal/HAWKI.Ks". The home cache is a valid read location for SVO
    # (that is where downloads land), so include it here.
    relpath = Path(config.svo_dir, *filter_name.split("/"))

    try:
        local_path = find_cached_file(
            relpath, config.package_name, extra_dirs=[__data_dir__]
        )
    except FileNotFoundError:
        target = config.write_cache_dir / relpath
        origin = (
            "http://svo2.cab.inta-csic.es/theory/fps3/"
            f"fps.php?ID={filter_name}"
        )
        local_path = Path(
            pooch.retrieve(
                origin, known_hash=None, fname=target.name, path=target.parent
            )
        )

    tbl = Table.read(local_path, format="votable")
    return tbl["Wavelength"], tbl["Transmission"]
