# -*- coding: utf-8 -*-
"""Utility functions for downloads. Subject to change."""

import warnings
from pathlib import Path

import pooch

from .configuration import config


__all__ = ["retriever",
           "download_file",
           "download_svo_filter"]


retriever = pooch.create(path=config.cache_dir,
                         base_url=config.database_url,
                         retry_if_failed=config.retry)
retriever.load_registry(config.registry_file)


def download_file(remote_url, local_name):
    """For backwards compatibility only."""
    warnings.warn("The download_file function is deprecated and will be "
                  "removed in v1.0. Please use retriever.fetch instead.",
                  DeprecationWarning, stacklevel=2)
    file = pooch.retrieve(remote_url, known_hash=None,
                          fname=local_name.name, path=local_name.parent)
    return file


def download_svo_filter(filter_name):
    """
    Query the SVO service for the true transmittance for a given filter.

    Parameters
    ----------
    filter_name : str
        Name of the filter as available on the spanish VO filter service
        e.g: ``Paranal/HAWKI.Ks``

    Returns
    -------
    filt_curve : tuple with wave and trans values
    """
    from astropy.table import Table

    data_dir = config.cache_dir

    origin = ("http://svo2.cab.inta-csic.es/theory/fps3/"
              f"fps.php?ID={filter_name}")

    local_path_cache = Path(data_dir, "svo_filters", filter_name)
    local_path_package = config.cache_dir / "svo" / filter_name

    if local_path_package.exists():
        local_path = local_path_package
    else:
        local_path = local_path_cache
        if not local_path.exists():
            download_file(origin, local_path)

    # raises ValueError if table is malformed
    # this can be used to catch problmes
    tbl = Table.read(local_path, format="votable")

    wave = tbl["Wavelength"]
    trans = tbl["Transmission"]

    return wave, trans
