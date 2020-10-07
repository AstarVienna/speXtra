# -*- coding: utf-8 -*-
"""
Database for speXtra
"""
import inspect
from posixpath import join as urljoin
import os
from urllib.error import URLError, HTTPError

import yaml
import tynt

from .utils import get_rootdir, database_url, download_file, dict_generator

# Configurations

__all__ = ["Database", "SpecLibrary", "SpectrumContainer", "ExtCurvesLibrary",
           "ExtCurveContainer", "FilterSystem", "FilterContainer",
           "get_filter_names", "get_filter_systems"]

# Do we need this?
__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# Default filters
with open(os.path.join(__data_dir__,  "default_filters.yml")) as filter_file:
    FILTER_DEFAULTS = yaml.safe_load(filter_file)

# Tables are displayed with a jsviewer by default
# Table.show_in_browser.__defaults__ = (5000, True, 'default', {'use_local_files': True},
#                                              None, 'display compact', None, 'idx')


class DataContainer:
    """
    Just a general container for the data

    It reads a YAML file and creates attributes with the kyes of the dictionary.

    Additionally it contains methods for nicely displaying the contents on the screen.
    """
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename) as f:
            self.meta = yaml.safe_load(f)
        for key in self.meta:
            setattr(self, key, self.meta[key])

    def dump(self):
        """
        (Try to) nicely dump the contents of the library

        Returns
        -------

        """
        print(yaml.dump(self.meta,
                        indent=4, sort_keys=False, default_flow_style=False))


class Database(DataContainer):
    """
    This class contains the database.


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

    def __init__(self, rootdir=get_rootdir(), remote_root=database_url()):

        if not remote_root.endswith('/'):
            remote_root = remote_root + '/'

        self.rootdir = rootdir
        self.remote_root = remote_root
        self.filename = "index.yml"
        self.path = self.abspath(self.filename)

        super().__init__(filename=self.path)

    def abspath(self, relpath):
        """
        Return absolute path to file or directory, ensuring that it exists.
        If it doesn't exist it will download it from the remote host

        Otherwise, just look for ``{relpath}``.
        """

        abspath = os.path.join(self.rootdir, relpath)
        if os.path.exists(abspath) is False:
            url = urljoin(self.remote_root, relpath)
            download_file(url, abspath)

        return abspath

    @property
    def pathlist(self):
        a = list(dict_generator(self.meta))
        separator = '/'

        return [separator.join(e) for e in a]


class SpecLibrary(DataContainer):
    """
    This class contains the information of a library

    """

    def __init__(self, library_name):

        database = Database()
        if library_name not in database.libraries:
            raise ValueError("library '%s' is not in the database" % library_name)

        self.directory = os.path.join("libraries", library_name)
        self.filename = "index.yml"
        self.relpath = os.path.join(self.directory, self.filename)
        self.path = database.abspath(self.relpath)

        super().__init__(filename=self.path)

        self.template_names = list(self.templates.keys())
        self.template_comments = list(self.templates.values())

    def __repr__(self):
        description = "Spectral Library: " + self.library_name + " " + self.title
        spec_cov = "spectral coverage: " + str(self.spectral_coverage)
        units = "wave_unit: " + self.wave_unit + "  flux_unit: " + self.flux_unit
        templates = "Templates: " + str(self.template_names)

        return ' %s \n %s \n %s \n %s' % (description, spec_cov, units, templates)


class FilterSystem(DataContainer):
    """
    This class contains the information of a filter system

    TODO: Quit tynt dependency
    """
    def __init__(self, filter_system):
        database = Database()
        if filter_system not in database.filter_systems:
            raise ValueError("filter system '%s' is not in the database" % filter_system)

        self.directory = os.path.join("filter_systems", filter_system)
        self.filename = "index.yml"
        self.relpath = os.path.join(self.directory, self.filename)
        self.path = database.abspath(self.relpath)

        super().__init__(filename=self.path)

        self.filter_names = list(self.filters.keys())
        self.filter_comments = list(self.filters.values())

    def __repr__(self):
        pass


class ExtCurvesLibrary(DataContainer):
    """
    Class that contains the information of the a Extinction Curve Library

    """
    def __init__(self, extinction_curve):
        database = Database()
        if extinction_curve not in database.extinction_curves:
            raise ValueError("extinction curve '%s' is not in the database" % extinction_curve)

        self.directory = os.path.join("extinction_curves", extinction_curve)
        self.filename = "index.yml"
        self.relpath = os.path.join(self.directory, self.filename)
        self.path = database.abspath(self.relpath)

        super().__init__(filename=self.path)

        self.curve_names = list(self.curves.keys())
        self.curve_comments = list(self.curves.values())

    def __repr__(self):
        pass


class SpectrumContainer(SpecLibrary):
    """
    A holder of spectral template information, including template characteristics and location
    """

    def __init__(self, template):
        self.template = template
        self.template_name = self.template.split("/").pop()

        library_name = self.template.rstrip("/" + self.template_name)
        super().__init__(library_name=library_name)

        if self.template_name not in self.template_names:
            raise ValueError("Template '%s' not in library" % self.template_name)

        self.template_comment = self.templates[self.template_name]
        self.filename = self.template_name + self.file_extension
        self.relpath = os.path.join(self.directory, self.filename)
        self.path = self.get_path()

        self._update_attrs()

    def get_path(self):
        database = Database()
        return database.abspath(self.relpath)

    def _update_attrs(self):
        """
        Here to just delete unnecessary stuff
        Returns
        -------
        """
        self.meta.pop("templates", None)
        self.meta.pop("summary", None)
        self.__delattr__("templates")
        self.__delattr__("summary")
        self.__delattr__("template_names")
        self.__delattr__("template_comments")

    def __repr__(self):
        s1 = "Spectral template: " + self.template_name
        return s1


class FilterContainer(FilterSystem):

    def __init__(self, filter):

        self.filter = filter
        if self.filter in FILTER_DEFAULTS:
            self.filter = FILTER_DEFAULTS[filter]

        if "/" not in self.filter_name:
            raise ValueError("not a valid filter %s" % self.filter_name)

        self.filter_name = self.filter.split("/").pop()
        filter_system = self.filter.rstrip("/" + self.filter_name)

        super().__init__(filter_system=filter_system)

        if self.filter_name not in self.filters:
            raise ValueError("Filter '%s' not in library", self.filter_name)

        self.filter_comment = self.filters[self.filter_name]
        self.filename = self.filter_name + self.file_extension
        self.relpath = os.path.join(self.directory, self.filename)
        self.path = self.get_path()

        self._update_attrs()

    def get_path(self):
        database = Database()
        return database.abspath(self.relpath)

    def _update_attrs(self):
        """
        Here to just delete unnecessary stuff
        """
        self.meta.pop("filterrs", None)
        self.meta.pop("summary", None)
        self.__delattr__("filters")
        self.__delattr__("summary")
        self.__delattr__("filters_names")
        self.__delattr__("filters_comments")

    #        download_file('http://svo2.cab.inta-csic.es/'
    #                      'theory/fps3/fps.php?ID={}'.format(self.filter_name),
    #                      os.path.join(database.rootdir, relpath))
    #        path = database.abspath(relpath)

    def __repr__(self):
        s = "Filter: " + self.filter_name
        return s


class ExtCurveContainer(ExtCurvesLibrary):

    def __init__(self, curve):
        self.curve = curve
        self.curve_name = self.curve.split("/").pop()

        name = self.template.rstrip("/" + self.curve_name)
        super().__init__(extinction_curve=name)

        if self.curve_name not in self.curve_names:
            raise ValueError("Extinction Curve '%s' not in library" % self.curve_name)

        self.curve_comment = self.curve_names[self.curve_name]
        self.filename = self.curve_name + self.file_extension
        self.relpath = os.path.join(self.directory, self.filename)
        self.path = self.get_path()

        self._update_attrs()

    def get_path(self):
        database = Database()
        return database.abspath(self.relpath)

    def _update_attrs(self):
        """
        Here to just delete unnecessary stuff
        Returns
        -------
        """
        self.meta.pop("curves", None)
        self.meta.pop("summary", None)
        self.__delattr__("curves")
        self.__delattr__("summary")
        self.__delattr__("curve_names")
        self.__delattr__("curve_comments")

    def __repr__(self):
        s = "Extinction curve: " + self.curve_name
        return s

#


def database_as_table():
    """
    Show the contents of the database as table
    Returns
    -------

    """

    pass


def database_as_tree():
    """
    Show the contents od the database as a tree
    Returns
    -------

    """
    pass


#    This is based on scopesim.effects.ter_curves_utils.py


def get_filter_systems():
    """
    Return a set of the different filter system available

    Returns
    -------

    """
    filters = tynt.FilterGenerator().available_filters()
    systems = {f.split("/")[0] for f in filters}
    return systems


def get_filter_names(system=None):
    """
    This function just returns the filters available from tynt
    if system= None returns all

    Returns
    -------

    """
    filter_list = tynt.FilterGenerator().available_filters()
    ord_list = [[f for f in filter_list if s in f] for s in get_filter_systems()]
    flat_list = [item for sublist in ord_list for item in sublist]

    if system is not None:
        flat_list = [f for f in filter_list if system in f]

    return flat_list
