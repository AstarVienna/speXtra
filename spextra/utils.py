# -*- coding: utf-8 -*-
"""Utility functions for SpeXtra."""

import numpy as np
import astropy.units as u



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

