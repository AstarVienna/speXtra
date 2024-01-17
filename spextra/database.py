# -*- coding: utf-8 -*-
"""Database for speXtra."""

import shutil
from io import StringIO
from pathlib import Path
from dataclasses import dataclass

import yaml

from astar_utils import NestedMapping

from .downloads import retriever
from .configuration import config

__all__ = ["load_yamldict", "spextra_database", "DEFAULT_DATA"]


class Database(NestedMapping):
    """Contains the database."""

    def __init__(self):
        self.database_url = config.database_url
        self.data_dir = config.cache_dir

        # TODO: change this to yaml constructor once in astar_utils
        super().__init__(load_yamldict("index.yml"))

    def fetch(self, filename: str) -> Path:
        """
        Return absolute path to file or directory, ensuring that it exists.

        If it doesn't exist it will download it from the remote host.
        Otherwise, just look for `relpath`.
        """
        filename = str(filename)
        abspath = self.data_dir / filename

        # don't download folder (.is_dir doesn't work here)
        if not abspath.suffix:
            return abspath

        abspath = Path(retriever.fetch(filename, progressbar=True))
        return abspath

    def clear_cache(self) -> None:
        """Remove local copies of all files in the database."""
        try:
            shutil.rmtree(self.data_dir)
            print(f"database at {self.data_dir} removed")
        except FileNotFoundError:
            print(f"database at {self.data_dir} doesn't exist")
            raise

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"

    def __str__(self) -> str:
        with StringIO() as str_stream:
            str_stream.write("Spextra Database:\n"
                             f"  Remote URL: {self.database_url}\n"
                             f"  Local path: {self.data_dir}\n")
            self.write_string(str_stream)
            output = str_stream.getvalue()
        return output


@dataclass
class DefaultData:
    """Small class just to define defaults spectra and filters."""

    filters: dict
    extcurves: dict
    spectra: dict

    @classmethod
    def from_yamls(cls):
        """Construct instance from default yaml files."""
        yamlfiles = ["default_filters.yml",
                     "default_spectra.yml",
                     "default_curves.yml"]
        paths = (Path(retriever.fetch(fname)) for fname in yamlfiles)
        yamls = dict_from_yamls(*paths)

        return cls(yamls["default_filters"],
                   yamls["default_curves"],
                   yamls["default_spectra"])


def load_yamldict(filename: str) -> dict:
    """Fetch YAML file from database and load contents into dict."""
    filepath = retriever.fetch(filename, progressbar=True)
    with open(filepath, encoding="utf-8") as file:
        yamldict = yaml.safe_load(file)
    return yamldict


def dict_from_yamls(*yamls):
    """Construct a nested dict from one or more yaml files."""
    outdict = {}
    for path in yamls:
        with path.open("r", encoding="utf-8") as file:
            outdict[path.stem] = yaml.safe_load(file)
    return outdict


DEFAULT_DATA = DefaultData.from_yamls()
spextra_database = Database()
