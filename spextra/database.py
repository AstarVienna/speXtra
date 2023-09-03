# -*- coding: utf-8 -*-
"""Database for speXtra."""

import shutil
from posixpath import join as urljoin
from pathlib import Path

import yaml

from .utils import Config, download_file, dict_generator

# Configurations

__all__ = ["Database",
           "SpecLibrary",
           "SpectrumContainer",
           "ExtCurvesLibrary",
           "ExtCurveContainer",
           "FilterSystem",
           "FilterContainer",
           "DefaultData"]

# Tables are displayed with a jsviewer by default
# Table.show_in_browser.__defaults__ = (5000, True, 'default', {'use_local_files': True},
#                                              None, 'display compact', None, 'idx')

CONF = Config()


class DataContainer:
    """
    Just a general container for the data.

    Reads a YAML file and creates attributes with the kyes of the dictionary.

    Additionally contains methods for nicely displaying the contents on screen.
    """
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename, encoding="utf-8") as file:
            self.meta = yaml.safe_load(file)
        for key in self.meta:
            setattr(self, key, self.meta[key])

    def dump(self) -> None:
        """
        (Try to) nicely dump the contents of the library

        Returns
        -------

        """
        print(yaml.dump(
            self.meta, indent=4, sort_keys=False, default_flow_style=False))


class Database(DataContainer):
    """
    Contains the database.

    It also acts as a lazy fetcher for remote data.
    When asked for local absolute path to a file or directory, Database
    checks if the file or directory exists locally and, if so, returns it.
    If it doesn't exist, it first determines where to get it from.
    It first downloads the file ``{remote_root}/redirects.json`` and checks
    it for a redirect from ``{relative_path}`` to a full URL. If no redirect
    exists, it uses ``{remote_root}/{relative_path}`` as the URL.
    It downloads then downloads the URL to ``{rootdir}/{relative_path}``.
    For directories, ``.tar.gz`` is appended to the
    ``{relative_path}`` before the above is done and then the
    directory is unpacked locally.
    Parameters
    ----------
    rootdir : str or callable
        The local root directory, or a callable that returns the local root
        directory given no parameters. (The result of the call is cached.)
        Using a callable allows one to customize the discovery of the root
        directory (e.g., from a config file), and to defer that discovery
        until it is needed.
    remote_root : str
        Root URL of the remote server.
    """

    def __init__(self, silent=False):
        self.database_url = CONF.get_database_url()
        self.data_dir = CONF.get_data_dir()
        self.ymlfile = "index.yml"
        self.path = self.abspath(self.ymlfile, silent=silent)

        super().__init__(filename=self.path)

        self.liblist, self.relpathlist = self._makelists()

    def abspath(self, relpath, reload=False, silent=False):
        """
        Return absolute path to file or directory, ensuring that it exists.

        If it doesn't exist it will download it from the remote host.
        Otherwise, just look for `relpath`.
        """
        # TODO: The name of this method does not at all suggest it will
        #       perform a download, which is non-trivial.
        relpath = Path(relpath)
        abspath = self.data_dir / relpath
        if reload or not abspath.exists():
            print(f"updating/loading '{relpath!s}'" )
            url = urljoin(self.database_url, *relpath.parts)
            download_file(url, str(abspath), silent=silent)

        return abspath

    def remove_database(self) -> None:
        """Remove database and all files."""
        try:
            shutil.rmtree(self.data_dir)
            print(f"database at {self.data_dir} removed")
        except FileNotFoundError:
            print(f"database at {self.data_dir} doesn't exist")
            raise

    def update(self) -> None:
        """Force an update of the cached database file."""
        self.abspath(self.ymlfile, reload=True,  silent=False)

    def _makelists(self):
        """Make lists of the libraries and paths in the database."""
        meta_list = list(dict_generator(self.meta))
        separator = "/"
        libs = [e[1:] for e in meta_list]
        liblist = [separator.join(e) for e in libs]
        relpathlist = [separator.join(e) for e in meta_list]

        return liblist, relpathlist

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"

    def __str__(self) -> str:
        outstr = ("Database:\n"
                  f"url: {self.database_url}\n"
                  f"path: {self.data_dir}\n")
        return outstr


class DefaultData:
    """Small class just to define defaults spectra and filters."""
    # TODO: where is this class used?? is it supposed to be inherited or smth?

    def __init__(self, **kwargs):

        database = Database()
        self.default_filters_file = database.abspath("default_filters.yml")
        self.default_spectra_file = database.abspath("default_spectra.yml")
        self.default_curves_file = database.abspath("default_curves.yml")

        self.filters = self.setdict(self.default_filters_file, **kwargs)
        self.spectra = self.setdict(self.default_spectra_file, **kwargs)
        self.extcurves = self.setdict(self.default_curves_file, **kwargs)

    @staticmethod
    def setdict(ymlfile, **kwargs):
        database = Database()
        path = database.abspath(ymlfile, **kwargs)

        with open(path, encoding="utf-8") as filename:
            dictionary = yaml.safe_load(filename)
        return dictionary

    @staticmethod
    def update() -> None:
        DefaultData(silent=False, reload=True)


class Library(DataContainer):
    """
    Class that contains the information of a Library.
    
    Either spectral of a filter system, extinction curves, etc.
    """
    def __init__(self, library_name):
        self.name = library_name
        self.title = "<untitled>"
        self.wave_unit = None
        self.file_extension = None
        self.items = []

        database = Database()
        if library_name not in database.liblist:
            raise ValueError(f"Library '{library_name}' not in the database")

        index = database.liblist.index(library_name)
        self.relpath = database.relpathlist[index]

        self.ymlfile = "index.yml"
        self.path = database.abspath(Path(self.relpath, self.ymlfile))
        self.dir = database.abspath(self.relpath)
        self.url = urljoin(database.data_dir, self.relpath, self.ymlfile)
        super().__init__(filename=self.path)

    def download_all(self) -> None:
        """Download the whole library."""
        database = Database()
        for item in self.items:
            database.abspath(Path(self.name, item))

    def remove(self) -> None:
        """Remove library and all files."""
        try:
            shutil.rmtree(self.dir)
            print(f"library {self.name} removed")
        except FileNotFoundError:
            print(f"library {self.name} doesn't exist")
            raise

    def update(self) -> None:
        """Update the library file."""
        database = Database()
        database.abspath(Path(self.relpath, self.ymlfile), reload=True,
                         silent=False)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.name!r})"


class SpecLibrary(Library):

    def __init__(self, library_name):
        self.templates = {}
        self.flux_unit = "<unknown>"
        self.spectral_coverage = []
        super().__init__(library_name)

        self.template_names = list(self.templates.keys())
        self.template_comments = list(self.templates.values())
        self.files = [t + self.file_extension for t in self.template_names]

        self.items = self.template_names

    def __str__(self) -> str:
        outstr = (f"Spectral Library '{self.name}': {self.title}\n"
                  f"spectral coverage: {', '.join(self.spectral_coverage)}\n"
                  f"wave_unit: {self.wave_unit}\n"
                  f"flux_unit: {self.flux_unit}\n"
                  f"Templates: {', '.join(self.template_names)}")
        return outstr


class FilterSystem(Library):
    """Contains the information of a filter system."""

    def __init__(self, filter_system):
        self.filters = {}
        self.spectral_coverage = []
        super().__init__(filter_system)

        self.filter_names = list(self.filters.keys())
        self.filter_comments = list(self.filters.values())
        self.files = [f + self.file_extension for f in self.filter_names]

        self.items = self.filter_names

    def __str__(self) -> str:
        filters = [f"{self.name}/{key}" for key in self.filters]
        outstr = (f"Filter system '{self.name}': {self.title}\n"
                  f"spectral coverage: {', '.join(self.spectral_coverage)}\n"
                  f"wave_unit: {self.wave_unit}\n"
                  f"filters: {', '.join(filters)}")
        return outstr


class ExtCurvesLibrary(Library):
    """Contains the information of the a Extinction Curve Library."""

    def __init__(self, curve_library):
        self.curves = {}
        super().__init__(curve_library)

        self.curve_names = list(self.curves.keys())
        self.curve_comments = list(self.curves.values())
        self.files = [e + self.file_extension for e in self.curve_names]

        self.items = self.curve_names

    def __str__(self) -> str:
        outstr = (f"Extinction Curves '{self.name}': {self.title}\n"
                  f"wave_unit: {self.wave_unit}\n"
                  f"Templates: {self.curves!s}")
        return outstr


class DBFile:
    """Mixin base class for database files."""

    _subclass_key = None

    def __init__(self, path=None):
        self.path = path
        self.meta = {}

    def remove(self) -> None:
        """Remove the file."""
        try:
            shutil.rmtree(self.path)
            print(f"library {self.path} removed")
        except FileNotFoundError:
            print(f"file {self.path} doesn't exist")
            raise

    def _update_attrs(self) -> None:
        """Delete unnecessary stuff."""
        self.meta.pop("summary", None)
        self.__dict__.pop("summary", None)
        self.__dict__.pop("ymlfile", None)

        if self._subclass_key is not None:
            self.meta.pop(f"{self._subclass_key}s", None)
            self.__dict__.pop(f"{self._subclass_key}s", None)
            self.__dict__.pop(f"{self._subclass_key}_names", None)
            self.__dict__.pop(f"{self._subclass_key}_comments", None)


class SpectrumContainer(SpecLibrary, DBFile):
    """
    Container of spectral template information.
    
    Including template characteristics and location.
    """

    _subclass_key = "template"

    def __init__(self, template):
        self.template = template
        *library_name, self.template_name = self.template.split("/")
        if isinstance(library_name, list):
            library_name = "/".join(library_name)
        super().__init__(library_name=library_name)

        if self.template_name not in self.template_names:
            raise ValueError(f"Template '{self.template}' not in library")

        self.datafile = self.template_name + self.file_extension

        database = Database()

        self.path = database.abspath(Path(self.relpath, self.datafile))
        self.url = urljoin(database.database_url, self.relpath, self.datafile)
        self.template_comment = self.templates[self.template_name]
        self.filename = self.path
        self._update_attrs()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.template!r})"

    def __str__(self) -> str:
        return f"Spectral template '{self.template}'"


class FilterContainer(FilterSystem, DBFile):

    _subclass_key = "filter"

    def __init__(self, filter_name):
        self.filter_name = filter_name
        *filter_system, self.basename = filter_name.split("/")
        if isinstance(filter_system, list):
            filter_system = "/".join(filter_system)

        super().__init__(filter_system=filter_system)
        if self.basename not in self.filters:
            raise ValueError(f"Filter {self.filter_name} not in library")

        database = Database()

        self.datafile = self.basename + self.file_extension
        self.path = database.abspath(Path(self.relpath, self.datafile))
        self.url = urljoin(database.database_url, self.relpath, self.datafile)
        self.filter_comment = self.filters[self.basename]
        self.filename = self.path

        self._update_attrs()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.filter_name!r})"

    def __str__(self) -> str:
        return f"Filter '{self.filter_name}'"


class ExtCurveContainer(ExtCurvesLibrary, DBFile):

    _subclass_key = "curve"

    def __init__(self, curve_name):
        self.curve_name = curve_name
        *curve_library, self.basename = curve_name.split("/")
        if isinstance(curve_library, list):
            curve_library = "/".join(curve_library)

        super().__init__(curve_library=curve_library)

        if self.basename not in self.curve_names:
            raise ValueError(f"Extinction Curve '{self.curve_name}' not in "
                             "library")

        database = Database()

        self.datafile = self.basename + self.file_extension
        self.path = database.abspath(Path(self.relpath, self.datafile))
        self.url = urljoin(database.database_url, self.relpath, self.datafile)
        self.curve_comment = self.curves[self.basename]
        self.filename = self.path

        self._update_attrs()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.curve_name!r})"

    def __str__(self) -> str:
        return f"Extinction curve '{self.curve_name}'"

# TODO: functions


def get_library(library_name):
    """Download library and all files."""
    raise NotImplementedError()


def get_filter_system(filter_system):
    """Download all filters of a system."""
    raise NotImplementedError()


def get_extcurves(extinction_curves):
    """Download all extinction curves."""
    raise NotImplementedError()


def download_database():
    """Download."""
    raise NotImplementedError()


def set_root_dir(path):
    """
    Set the root dir where the files will be stored.

    Probably best way is to write a small file in the data directory.

    Also get root dir should be also able to read from that file.

    Database should be get rid of initialiting parameters or use them
    to set new directories
    """
    raise NotImplementedError()


def database_as_table():
    """Show the contents of the database as table."""
    raise NotImplementedError()


def database_as_tree():
    """Show the contents od the database as a tree."""
    raise NotImplementedError()
