# -*- coding: utf-8 -*-
"""Contains database file class and container classes derived from it."""

import re
import warnings
from pathlib import Path

from .exceptions import NotInLibraryError, ConstructorError
from .database import spextra_database
from .libraries import Library, SpecLibrary, FilterSystem, ExtCurvesLibrary


__all__ = ["SpectrumContainer", "ExtCurveContainer", "FilterContainer"]

NAME_REGEX = re.compile(r"^(?P<libname>\w+(/\w+)*)/(?P<basename>[\w\.\-]+)$")


class DBFile:
    """Base class for database files."""

    library = None  # necessary to avoid access before init ran (REALLY??)
    _subclass_str = "DBFile"
    _subclass_library = Library

    def __init__(self, name: str):
        # To avoid name conflict with backwards-compatible meaning
        # In the future, just call it .name
        self.name = name
        name_match = NAME_REGEX.match(name)
        if not name_match:
            raise ConstructorError(
                f"{self._subclass_str} name '{name}' must match the pattern "
                "'library_name/file_name'.")

        self.basename = name_match["basename"]
        self.library = self._subclass_library(name_match["libname"])
        # TODO: maybe an info log if not in DB?

    @property
    def is_in_library(self) -> bool:
        """Return True if file is part of the parent library."""
        return self.basename in self.library

    @property
    def datafile(self):
        raise AttributeError(
            "The .datafile attribute is deprecated and wil be fully removed in"
            "version 0.45.0. Please use .filename instead.")

    @property
    def path(self) -> Path:
        """Path to the cached file."""
        if not self.is_in_library:
            raise NotInLibraryError(f"{self._subclass_str} '{self.basename}' "
                                    "not in library")
        return spextra_database.fetch(f"{self.library.path}/{self.filename}")

    @property
    def filename(self):
        """Name and extension of the file.

        TODO: versionchanged placeholder here!!
        """
        return self.basename + self.library.file_extension

    @property
    def description(self) -> str:
        """Description or comment of the file as defined in the library."""
        return self.library[self.basename]

    def remove(self) -> None:
        """Remove the file."""
        self.path.unlink()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.name!r})"

    def __str__(self) -> str:
        return f"{self._subclass_str} '{self.name}' ({self.description})"


class SpectrumContainer(DBFile):
    """Container for a spectral template file."""

    _subclass_str = "Spectral template"
    _subclass_library = SpecLibrary

    @property
    def template(self) -> str:
        """Name of the spectrum template."""
        return self.name

    @property
    def template_name(self):
        """Deprecated feature."""
        warnings.warn("The .template_name property is deprecated and will be "
                      "removed in v1.0. Please use the identical .basename "
                      "instead!", FutureWarning, 2)
        return self.basename

    @property
    def template_comment(self):
        """Deprecated feature."""
        warnings.warn("The .template_comment property is deprecated and will be "
                      "removed in v1.0. Please use the identical .description "
                      "instead!", FutureWarning, 2)
        return self.description


class FilterContainer(DBFile):
    """Container for a filter curve file."""

    _subclass_str = "Filter"
    _subclass_library = FilterSystem

    @property
    def filter_name(self) -> str:
        """Name of the filter."""
        return self.name

    @property
    def filter_comment(self):
        """Deprecated feature."""
        warnings.warn("The .filter_comment property is deprecated and will be "
                      "removed in v1.0. Please use the identical .description "
                      "instead!", FutureWarning, 2)
        return self.description


class ExtCurveContainer(DBFile):
    """Container for extinction curve file."""

    _subclass_str = "Extinction Curve"
    _subclass_library = ExtCurvesLibrary

    @property
    def curve_name(self) -> str:
        """Name of the extinction curve."""
        return self.name

    @property
    def curve_comment(self):
        """Deprecated feature."""
        warnings.warn("The .curve_comment property is deprecated and will be "
                      "removed in v1.0. Please use the identical .description "
                      "instead!", FutureWarning, 2)
        return self.description
