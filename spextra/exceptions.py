# -*- coding: utf-8 -*-
"""Spextra-specfic exceptions and errors, for encapsulation."""


__all__ = [
    "SpextraError",
    "NotInDatabaseError",
    "NotInLibraryError",
    "ArgumentError",
    "ConstructorError",
    ]


class SpextraError(Exception):
    """Base class for errors in this Package."""


class NotInDatabaseError(SpextraError, ValueError):
    """Requested library is not in the database."""


class NotInLibraryError(SpextraError, ValueError):
    """Requested file is not in the library."""


class ArgumentError(SpextraError, ValueError):
    """Error in passed arguments."""


class ConstructorError(ArgumentError):
    """Error in constructor arguments."""


class UnitError(SpextraError):
    """Unit not understood by SYNPHOT or astropy."""
