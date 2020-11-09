# -*- coding: utf-8 -*-
"""
Database for speXtra
"""
import inspect
from posixpath import join as urljoin
import os
from urllib.error import URLError, HTTPError

import yaml

from .utils import get_rootdir, database_url, download_file, dict_generator

# Configurations

__all__ = ["Database", "SpecLibrary", "SpectrumContainer", "ExtCurvesLibrary",
           "ExtCurveContainer", "FilterSystem", "FilterContainer"]

# Do we need this?
__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# Default filters

# Tables are displayed with a jsviewer by default
# Table.show_in_browser.__defaults__ = (5000, True, 'default', {'use_local_files': True},
#                                              None, 'display compact', None, 'idx')

with open(os.path.join(__data_dir__,  "default_filters.yml")) as filter_file:
    FILTER_DEFAULTS = yaml.safe_load(filter_file)


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
        self.ymlfile = "index.yml"
        self.path = self.abspath(self.ymlfile)

        super().__init__(filename=self.path)

        self._makelists()

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

    def _makelists(self):
        a = list(dict_generator(self.meta))
        separator = '/'
        libs = [e[1:] for e in a]
        self.liblist = [separator.join(e) for e in libs]
        self.relpathlist = [separator.join(e) for e in a]


class Library(DataContainer):
    """
    Class that contains the information of a Library, either spectral
    of a filter system, extinction curves, etc
    """
    def __init__(self, library_name):

        database = Database()

        if library_name not in database.liblist:
            raise ValueError("Library '%s' not in the database" % library_name)
        else:
            index = database.liblist.index(library_name)
            self.relpath = database.relpathlist[index]

        self.ymlfile = "index.yml"
        self.path = database.abspath(os.path.join(self.relpath, self.ymlfile))
        self.dir = database.abspath(self.relpath)
        self.url = urljoin(database.remote_root, self.relpath, self.ymlfile)
        super().__init__(filename=self.path)


class SpecLibrary(Library):

    def __init__(self, library_name):

        super().__init__(library_name)

        self.template_names = list(self.templates.keys())
        self.template_comments = list(self.templates.values())
        self.files = [t + self.file_extension for t in self.template_names]

    def __repr__(self):
        description = "Spectral Library: " + self.library_name + " " + self.title
        spec_cov = "spectral coverage: " + str(self.spectral_coverage)
        units = "wave_unit: " + self.wave_unit + "  flux_unit: " + self.flux_unit
        templates = "Templates: " + str(self.template_names)

        return ' %s \n %s \n %s \n %s' % (description, spec_cov, units, templates)


class FilterSystem(Library):
    """
        This class contains the information of a filter system

        """

    def __init__(self, filter_system):

        super().__init__(filter_system)

        self.filter_names = list(self.filters.keys())
        self.filter_comments = list(self.filters.values())
        self.files = [f + self.file_extension for f in self.filter_names]


class ExtCurvesLibrary(Library):
    """
    Class that contains the information of the a Extinction Curve Library
    """
    def __init__(self, curve_library):

        super().__init__(curve_library)

        self.curve_names = list(self.curves.keys())
        self.curve_comments = list(self.curves.values())
        self.files = [e + self.file_extension for e in self.curve_names]


class SpectrumContainer(SpecLibrary):
    """
    A holder of spectral template information, including template characteristics and location
    """

    def __init__(self, template):
        self.template = template
        self.template_name = os.path.basename(self.template)

        library_name = os.path.split(self.template)[0]
        super().__init__(library_name=library_name)

        if self.template_name not in self.template_names:
            raise ValueError("Template '%s' not in library" % self.template_name)

        self.datafile = self.template_name + self.file_extension

        database = Database()

        self.path = database.abspath(os.path.join(self.relpath, self.datafile))
        self.url = urljoin(database.remote_root, self.relpath, self.datafile)
        self.template_comment = self.templates[self.template_name]
        self.filename = self.path
        self._update_attrs()

    def _update_attrs(self):
        """
        Here to just delete unnecessary stuff
        Returns
        -------
        """
        self.meta.pop("templates", None)
        self.meta.pop("summary", None)
        self.__dict__.pop("templates", None)
        self.__dict__.pop("summary", None)
        self.__dict__.pop("template_names", None)
        self.__dict__.pop("template_comments", None)
        self.__dict__.pop("ymlfile", None)

    def __repr__(self):
        s1 = "Spectral template: " + self.template_name
        return s1


class FilterContainer(FilterSystem):

    def __init__(self, filter_name):

        self.filter_name = filter_name
        if self.filter_name in FILTER_DEFAULTS:
            self.filter_name = FILTER_DEFAULTS[filter_name]

        self.name = os.path.basename(self.filter_name)
        filter_system = os.path.split(self.filter_name)[0]

        super().__init__(filter_system=filter_system)

        if self.name not in self.filters:
            raise ValueError("Filter '%s' not in library", self.filter_name)

        database = Database()

        self.datafile = self.name + self.file_extension
        self.path = database.abspath(os.path.join(self.relpath, self.datafile))
        self.url = urljoin(database.remote_root, self.relpath, self.datafile)
        self.filter_comment = self.filters[self.name]
        self.filename = self.path

        self._update_attrs()

    def _update_attrs(self):
        """
        Here to just delete unnecessary stuff
        """
        self.meta.pop("filterrs", None)
        self.meta.pop("summary", None)
        self.__dict__.pop("filters", None)
        self.__dict__.pop("summary", None)
        self.__dict__.pop("filter_names", None)
        self.__dict__.pop("filter_comments", None)
        self.__dict__.pop("ymlfile", None)

    #        download_file('http://svo2.cab.inta-csic.es/'
    #                      'theory/fps3/fps.php?ID={}'.format(self.filter_name),
    #                      os.path.join(database.rootdir, relpath))
    #        path = database.abspath(relpath)

    def __repr__(self):
        s = "Filter: " + self.filter_name
        return s


class ExtCurveContainer(ExtCurvesLibrary):

    def __init__(self, curve_name):

        self.curve_name = curve_name
        self.cname = os.path.basename(self.curve_name)
        curve_library = os.path.split(self.curve_name)[0]
        print(self.curve_name, self.cname, curve_library)
        super().__init__(curve_library=curve_library)

        if self.cname not in self.curve_names:
            raise ValueError("Extinction Curve '%s' not in library" % self.curve_name)



        database = Database()

        self.datafile = self.cname + self.file_extension
        self.path = database.abspath(os.path.join(self.relpath, self.datafile))
        self.url = urljoin(database.remote_root, self.relpath, self.datafile)
        self.curve_comment = self.curves[self.cname]
        self.filename = self.path

        self._update_attrs()

    def _update_attrs(self):
        """
        Here to just delete unnecessary stuff
        """
        self.meta.pop("curves", None)
        self.meta.pop("summary", None)
        self.__dict__.pop("curves", None)
        self.__dict__.pop("summary", None)
        self.__dict__.pop("curve_names", None)
        self.__dict__.pop("curve_comments", None)
        self.__dict__.pop("ymlfile", None)

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



