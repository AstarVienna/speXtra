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
    """Base class for database files.


    .. versionchanged:: 0.44.0

       Attempting to construct an instance that does not point to a valid file
       in the specified library will now immediately raise an error, rather
       than only when ``.path`` is accessed the first time. The exception to
       this are filters, which might be found in the SVO.
    """

    library = None  # necessary to avoid access before init ran (REALLY??)
    _subclass_str = "DBFile"
    _subclass_library = Library
    _strict = True  # overridded for Filters -> SVO fallback

    def __init__(self, name: str):
        self.name = name
        name_match = NAME_REGEX.match(name)
        if not name_match:
            raise ConstructorError(
                f"{self._subclass_str} name '{name}' must match the pattern "
                "'library_name/file_name'.")

        self.basename = name_match["basename"]
        self.library = self._subclass_library(name_match["libname"])

        if self._strict and not self.is_in_library:
            raise NotInLibraryError(
                f"{self._subclass_str} '{self.basename}' "
                f"not in library '{self.library.name}'."
            )

    @property
    def is_in_library(self) -> bool:
        """Return True if file is part of the parent library."""
        return self.basename in self.library

    @property
    def datafile(self):
        """
        .. deprecated:: 0.44.0
        """
        raise AttributeError(
            "The .datafile attribute is deprecated and wil be fully removed in"
            "version 0.45.0. Please use .filename instead.")

    @property
    def path(self) -> Path:
        """Path to the cached file."""
        if not self.is_in_library:
            raise NotInLibraryError(
                f"{self._subclass_str} '{self.basename}' "
                f"not in library '{self.library.name}'."
            )
        return spextra_database.fetch(f"{self.library.path}/{self.filename}")

    @property
    def filename(self):
        """Name and extension of the file.

        .. versionchanged:: 0.44.0

           This attribute used to return ``self.library.name`` before v0.44.0,
           it now returns ``self.basename + self.library.file_extension``.
           If you need the name of the underlying library, use
           ``.library.name`` directly.
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
        """Deprecated feature.

        .. deprecated:: 0.44.0

           Please use the identical ``.basename`` instead!
        """
        warnings.warn("The .template_name property is deprecated and will be "
                      "removed in v1.0. Please use the identical .basename "
                      "instead!", FutureWarning, 2)
        return self.basename

    @property
    def template_comment(self):
        """Deprecated feature.

        .. deprecated:: 0.44.0

           Please use the identical ``.description`` instead!
        """
        warnings.warn("The .template_comment property is deprecated and will be "
                      "removed in v1.0. Please use the identical .description "
                      "instead!", FutureWarning, 2)
        return self.description


class FilterContainer(DBFile):
    """Container for a filter curve file."""

    _subclass_str = "Filter"
    _subclass_library = FilterSystem
    _strict = False

    @property
    def filter_name(self) -> str:
        """Name of the filter."""
        return self.name

    @property
    def filter_comment(self):
        """Deprecated feature.

        .. deprecated:: 0.44.0

           Please use the identical ``.description`` instead!
        """
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
        """Deprecated feature.

        .. deprecated:: 0.44.0

           Please use the identical ``.description`` instead!
        """
        warnings.warn("The .curve_comment property is deprecated and will be "
                      "removed in v1.0. Please use the identical .description "
                      "instead!", FutureWarning, 2)
        return self.description
