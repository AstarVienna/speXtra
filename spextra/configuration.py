# -*- coding: utf-8 -*-
"""Holds the Config class and an instance to use as global configuration."""

from pathlib import Path
from dataclasses import dataclass, asdict
import yaml

from .exceptions import SpextraError

__all__ = ["config"]
__data_dir__ = Path(__file__).parent / "data"
__config_file__ = __data_dir__ / "config.yml"


@dataclass
class Config:
    """
    Set up the configuration for the database.

    This class reads a yaml file in the directory data_dir that contains the
    default configuration of the database that contain default values for the
    following attributes:

    database_url: The location of the remote database
    cache_dir: The location of the database in the local hard disk, where the
        files are saved. Default is Null (None) which then lead to
        .spextra_cache
    retry: The number of times to re-attempt a download.

    Normally, the user will not need to call this class, except when chaning

    If the user needs to change these values the class can be simply called as

    >> Config(cache_dir="path_to_new_data_directory")

    and the new path will be used from that moment on
    """

    database_url: str
    cache_dir: Path
    retry: int
    registry_file: Path

    def __post_init__(self):
        self.cache_dir = Path(self.cache_dir)
        self.registry_file = Path(self.registry_file)

        if not self.cache_dir.is_absolute():
            self.cache_dir = Path.home() / self.cache_dir

        if not self.registry_file.is_absolute():
            self.registry_file = __data_dir__ / self.registry_file

        if not self.registry_file.exists():
            raise SpextraError("registry file not found")

        # TODO: dump back into config file
        # TODO: add method to restore defaults

    @classmethod
    def from_yaml(cls, path: Path):
        with path.open(encoding="utf-8") as file:
            return cls(**yaml.safe_load(file))

    def __str__(self) -> str:
        sout = "SpeXtra Configuration:\n"
        for k, v in asdict(config).items():
            sout += f"  {k}: {v!s}\n"
        return sout


config = Config.from_yaml(__config_file__)
