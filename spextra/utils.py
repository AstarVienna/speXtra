# -*- coding: utf-8 -*-
"""Utility functions for SpeXtra."""

import astropy.units as u


def _angstrom_value(value):
    if isinstance(value, u.Quantity):
        return value.to(u.AA, equivalencies=u.spectral()).value
    return value
