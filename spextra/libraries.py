# -*- coding: utf-8 -*-
"""Contains ldatabase library and subclasses."""

import warnings
import shutil
from pathlib import Path
from collections.abc import Mapping

from synphot.units import validate_unit
from synphot.exceptions import SynphotError

from .exceptions import UnitError, ConstructorError
from .database import spextra_database, load_yamldict


__all__ = ["SpecLibrary", "ExtCurvesLibrary", "FilterSystem"]


class Library(Mapping):
    """
    Class that contains the information of a Library.

    Either spectral of a filter system, extinction curves, etc.

    Base class, not meant to be instantiated directly.
    """

    db_dir = None
    aliases = {
        "items": "items",
        "wave_col": "wave_col",
        "flux_col": "flux_col",
        "flux_unit": "flux_unit",
        }

    def __init__(self, name):
        self.name = name

        if self.is_in_database:
            yamldict = load_yamldict(f"{self.path}/index.yml")
        else:
            yamldict = {}

        self.title = yamldict.get("title", "<untitled>")
        self._items = yamldict.get(self.aliases["items"], {})
        self.spectral_coverage = yamldict.get("spectral_coverage", [])

        self.read_kwargs = self._get_read_kwargs(yamldict)
        self._validate_units()

        self.file_extension = yamldict.get("file_extension", None)

    def _get_read_kwargs(self, yamldict) -> dict:
        read_kwargs = {
            "wave_col": yamldict.get(self.aliases["wave_col"], "WAVELENGTH"),
            "flux_col": yamldict.get(self.aliases["flux_col"], "FLUX"),
            "wave_unit": yamldict.get("wave_unit", "Angstrom"),
            "flux_unit": yamldict.get(self.aliases["flux_unit"], "FLAM"),
        }
        return read_kwargs

    def _validate_units(self) -> None:
        for key in ("wave_unit", "flux_unit"):
            value = self.read_kwargs[key]
            # HACK to understand ext. units
            if value == "Av/E(B-V)":
                value = "mag"
            try:
                self.read_kwargs[key] = validate_unit(value)
            except (SynphotError, ValueError) as err:
                raise UnitError(f"{value} not understood.") from err

    @property
    def is_in_database(self) -> bool:
        """Return True if library is defined in the database instance."""
        return f"!{self.db_dir}.{self.name}" in spextra_database

    @property
    def path(self):
        return f"{self.db_dir}/{self.name}"

    @property
    def url(self):
        pass  # Not yet supported

    @property
    def files(self):
        """Get a list of all filenames in the library."""
        return [f"{fname}{self.file_extension}" for fname in self]

    @property
    def data_type(self) -> str:
        """Determine data type based on known file extensions."""
        if self.file_extension == ".fits":
            return "fits"
        if self.file_extension == ".dat":
            return "ascii"
        return self.file_extension

    def download_all(self) -> None:
        """Download the whole library."""
        # database = Database()
        for key in self:
            spextra_database.fetch(Path(self.name, key))

    def clear_cache(self) -> None:
        """Remove local copies of all files in the library."""
        # database = Database()
        lib_dir = spextra_database.data_dir / self.path
        try:
            shutil.rmtree(lib_dir)
            print(f"library {self.name} removed")
        except FileNotFoundError:
            print(f"library {self.name} doesn't exist")
            raise

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.name!r})"

    def __getitem__(self, key):
        return self._items[key]

    def __iter__(self):
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)


class SpecLibrary(Library):
    """Library of spectral templates."""

    db_dir = "libraries"
    aliases = {
        "items": "templates",
        "wave_col": "wave_column_name",
        "flux_col": "flux_column_name",
        "flux_unit": "flux_unit",
        }

    def __init__(self, name=None, library_name=None):
        if library_name is not None:
            warnings.warn("Constructing using the 'library_name' argument is "
                          "deprecated and will raise an error in v1.0. Please "
                          "use the more general 'name' argument.",
                          DeprecationWarning, 2)
            name = library_name
        if name is None and library_name is None:
            raise ConstructorError("name must be passed to constructor")

        super().__init__(name)

    @property
    def library_name(self):
        """Deprecated feature."""
        warnings.warn("The alias .library_name is deprecated and will be "
                      "removed in v1.0. Please use .name instead!",
                      DeprecationWarning, 2)
        return self.name

    @property
    def templates(self):
        """Deprecated feature."""
        warnings.warn("The alias .templates is deprecated and will be "
                      "removed in v1.0. Please use .items() instead!",
                      DeprecationWarning, 2)
        return list(self.items())

    @property
    def template_names(self):
        """Deprecated feature."""
        warnings.warn("The alias .template_names is deprecated and will be "
                      "removed in v1.0. Please use the  more general "
                      ".keys() instead!", DeprecationWarning, 2)
        return list(self.keys())

    @property
    def template_comments(self):
        """Deprecated feature."""
        warnings.warn("The properts .template_comments is deprecated and will "
                      "be removed in v1.0.", DeprecationWarning, 2)
        return list(self.values())

    def __str__(self) -> str:
        outstr = (f"Spectral Library '{self.name}': {self.title}\n"
                  f"  spectral coverage: {', '.join(self.spectral_coverage)}\n"
                  f"  wave_unit: {self.read_kwargs['wave_unit']}\n"
                  f"  flux_unit: {self.read_kwargs['flux_unit']}\n"
                  f"  Templates: {', '.join(self.keys())}")
        return outstr


class FilterSystem(Library):
    """Library of filters comprising a filter system."""

    db_dir = "filter_systems"
    aliases = {
        "items": "filters",
        "wave_col": "wave_col",
        "flux_col": "flux_col",
        "flux_unit": "flux_unit",
        }

    def __init__(self, name=None, filter_system=None):
        if filter_system is not None:
            warnings.warn("Constructing using the 'filter_system' argument is "
                          "deprecated and will raise an error in v1.0. Please "
                          "use the more general 'name' argument.",
                          DeprecationWarning, 2)
            name = filter_system
        if name is None and filter_system is None:
            raise ConstructorError("name must be passed to constructor")

        super().__init__(name)

    @property
    def filter_system(self):
        """Deprecated feature."""
        warnings.warn("The alias .filter_system is deprecated and will be "
                      "removed in v1.0. Please use .name instead!",
                      DeprecationWarning, 2)
        return self.name

    @property
    def filters(self):
        """Deprecated feature."""
        warnings.warn("The alias .filters is deprecated and will be "
                      "removed in v1.0. Please use .items() instead!",
                      DeprecationWarning, 2)
        return list(self.items())

    @property
    def filter_names(self):
        """Deprecated feature."""
        warnings.warn("The alias .filter_names is deprecated and will be "
                      "removed in v1.0. Please use the  more general "
                      ".keys() instead!", DeprecationWarning, 2)
        return list(self.keys())

    @property
    def filter_comments(self):
        """Deprecated feature."""
        warnings.warn("The properts .filter_comments is deprecated and will "
                      "be removed in v1.0.", DeprecationWarning, 2)
        return list(self.values())

    def __str__(self) -> str:
        filters = [f"{self.name}/{key}" for key in self]
        outstr = (f"Filter system '{self.name}': {self.title}\n"
                  f"  spectral coverage: {', '.join(self.spectral_coverage)}\n"
                  f"  wave_unit: {self.read_kwargs['wave_unit']}\n"
                  f"  filters: {', '.join(filters)}")
        return outstr


class ExtCurvesLibrary(Library):
    """Library of extinction curves."""

    db_dir = "extinction_curves"
    aliases = {
        "items": "curves",
        "wave_col": "wave_column",
        "flux_col": "extinction_column",
        "flux_unit": "extinction_unit",
        }

    def __init__(self, name=None, curve_library=None):
        if curve_library is not None:
            warnings.warn("Constructing using the 'curve_library' argument is "
                          "deprecated and will raise an error in v1.0. Please "
                          "use the more general 'name' argument.",
                          DeprecationWarning, 2)
            name = curve_library
        if name is None and curve_library is None:
            raise ConstructorError("name must be passed to constructor")

        super().__init__(name)

    @property
    def curves(self):
        """Deprecated feature."""
        warnings.warn("The alias .curves is deprecated and will be "
                      "removed in v1.0. Please use .items() instead!",
                      DeprecationWarning, 2)
        return list(self.items())

    @property
    def curve_names(self):
        """Deprecated feature."""
        warnings.warn("The alias .curve_names is deprecated and will be "
                      "removed in v1.0. Please use the  more general "
                      ".keys() instead!", DeprecationWarning, 2)
        return list(self.keys())

    @property
    def curve_comments(self):
        """Deprecated feature."""
        warnings.warn("The properts .curve_comments is deprecated and will be "
                      "removed in v1.0.", DeprecationWarning, 2)
        return list(self.values())

    def __str__(self) -> str:
        outstr = (f"Extinction Curves '{self.name}': {self.title}\n"
                  f"  wave_unit: {self.read_kwargs['wave_unit']}\n"
                  f"  Templates: {self._items!s}")
        return outstr
