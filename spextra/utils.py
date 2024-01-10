# -*- coding: utf-8 -*-
"""Utility functions for SpeXtra."""

from urllib.request import urlopen, Request
from urllib.error import URLError

from pathlib import Path

import numpy as np
import astropy.units as u


__all__ = ["is_url"]

ARRAY_LIKE = (np.ndarray, list, tuple)

__data_dir__ = Path(__file__).parent / "data"


def _ensure_list(value):
    if isinstance(value, ARRAY_LIKE):
        return list(value)
    return [value]


def _angstrom_qty(value):
    if isinstance(value, u.Quantity):
        return value.to(u.AA, equivalencies=u.spectral())
    return value * u.AA


def _angstrom_value(value):
    if isinstance(value, u.Quantity):
        return _angstrom_qty(value).value
    return value


def _abmag_qty(value):
    if not isinstance(value, u.Quantity):
        return value * u.ABmag
    return value


def is_url(url: str) -> bool:
    """Check if connection works."""
    try:
        request = Request(url)
        request.get_method = lambda: "HEAD"
        urlopen(request)
    except URLError:
        return False
    except ValueError:
        return False
    return True
