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

from .utils import get_rootdir, database_url, download_file

# Configurations

__all__ = ["SpecDatabase", "SpecLibrary", "SpectralTemplate", "ExtCurvesLibrary",
           "ExtinctionCurve", "FilterSystem", "Filter",
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


class SpecDatabase:
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

    def __init__(self, rootdir=get_rootdir(),
                 remote_root=database_url()):
        if not remote_root.endswith('/'):
            remote_root = remote_root + '/'

        self.rootdir = rootdir
        self.remote_root = remote_root

        self.contents = self.get_yaml_contents("index.yml")
        self.libraries = [lib for lib in self.contents["libraries"]]
        self.extinction_curves = [ext for ext in self.contents["extinction_curves"]]
        self.filter_systems = [filt for filt in self.contents["filter_systems"]]

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

    def get_yaml_contents(self, relpath):
        """
        read a yaml file from a relative url

        Parameters
        ----------
        relpath: The relative path to a file either local or in remote host

        Returns
        -------
        dict with the contents of the yaml file
        """
        filename = self.abspath(relpath)
        with open(filename) as f:
            data = yaml.safe_load(f)

        return data

    def __repr__(self):
        rootdir = " local data directory: " + self.rootdir
        database_url = "data base url: " + self.remote_root
        libs = "libraries:" + ' ' + str(self.libraries)
        exts = "extinction curves:" + ' ' + str(self.extinction_curves)
        filts = "filter systems:" + ' ' + str(self.filter_systems)

        return '%s \n %s \n %s \n %s \n %s \n %s' % (rootdir, database_url,
                                                     'Database contents:', libs, exts, filts)


class SpecLibrary:
    """
    This class contains the information of a library

    """

    def __init__(self, name):
        self.name = name
        database = SpecDatabase()
        if self.name not in database.libraries:
            raise ValueError("library %s unknown, please check" % self.name)

        relpath = urljoin("libraries", self.name, "index.yml")

        self.meta = database.get_yaml_contents(relpath)

        self.location = urljoin(database.remote_root, relpath)
        self.library_name = self.meta["library_name"]
        self.title = self.meta["title"]
        self.type = self.meta["type"]
        self.summary = self.meta["summary"]
        self.reference = self.meta["reference"]
        self.link = self.meta["link"]
        self.spectral_coverage = self.meta["spectral_coverage"]
        self.resolution = self.meta["resolution"]
        self.wave_unit = self.meta["wave_unit"]
        self.flux_unit = self.meta["flux_unit"]
        self.wave_column_name = self.meta["wave_column_name"]
        self.flux_column_name = self.meta["flux_column_name"]
        self.data_type = self.meta["data_type"]
        self.file_extension = self.meta["file_extension"]
        self.templates = list(self.meta["templates"].keys())
        self.template_comments = [self.meta["templates"][k] for k in self.templates]

    def dump(self):
        """
        Nicely dump the contents of the library

        Returns
        -------

        """
        print(yaml.dump(self.meta,
                        indent=4, sort_keys=False, default_flow_style=False))

    def __repr__(self):
        description = "Spectral Library: " + self.name + " " + self.title
        spec_cov = "spectral coverage: " + str(self.spectral_coverage)
        units = "wave_unit: " + self.wave_unit + "  flux_unit: " + self.flux_unit
        templates = "Templates: " + str(self.templates)

        return ' %s \n %s \n %s \n %s' % (description, spec_cov, units, templates)


class SpectralTemplate:
    """
    A holder of spectral template information, including template characteristics and location
    """

    def __init__(self, template):
        self.library_name, self.template_name = template.split("/")
        library = SpecLibrary(self.library_name)
        if self.template_name not in library.templates:
            raise ValueError("Template not in library", self.template_name)
        self.resolution = library.resolution
        self.wave_unit = library.wave_unit
        self.flux_unit = library.flux_unit
        self.wave_column_name = library.wave_column_name
        self.flux_column_name = library.flux_column_name
        self.data_type = library.data_type
        self.file_extension = library.file_extension
        self.filename = self.template_name + self.file_extension
        self.path = self.get_path()

    def get_path(self):
        database = SpecDatabase()
        relpath = urljoin("libraries", self.library_name, self.filename)
        return database.abspath(relpath)

    @property
    def meta(self):
        library = SpecLibrary(self.library_name)
        del library.meta["templates"]
        del library.meta["summary"]
        return library.meta

    def __repr__(self):
        s1 = "Spectral template: " + self.template_name
        return s1


class FilterSystem:
    """
    ``FilterSystem`` holds all the information of a particular filter system which is in
    turn passed to ``Filter``

    Because we are currently using filters served by the SVO the attributes
    of these filter systems are not known and thus initialized with none.

    This might change in the future
    """

    def __init__(self, filter_system):
        self.filter_system = filter_system
        self.instrument = None
        self.title = None
        self.author = None
        self.source = None
        self.spectral_coverage = None
        self.wave_unit = None
        self.data_type = None
        self.file_extension = None
        self.filters = None
        self.meta = None

        relpath = urljoin("filter_systems", self.filter_system, "index.yml")
        database = SpecDatabase()
        if self.filter_system in database.filter_systems:
            self.meta = database.get_yaml_contents(relpath)

#        elif self.filter_system not in get_filter_systems():  # get_filter_systems does not report
#            data_dict = None  # Filter is in the SVO          # all filters at SVO only those created with tynt
#                                                              # if filter is not at SVO no error is raised
#        else:
#            raise ValueError("Filter system %s is unknown" % self.filter_system)

        self._update_atributes()

    def _update_atributes(self):

        if self.meta is not None:
            self.filter_system = self.meta["filter_system"]
            self.instrument = self.meta["instrument"]
            self.title = self.meta["title"]
            self.author = self.meta["author"]
            self.source = self.meta["source"]
            self.spectral_coverage = self.meta["spectral_coverage"]
            self.wave_unit = self.meta["wave_unit"]
            self.file_extension = self.meta["file_extension"]
            self.filters = list(self.meta["filters"].keys())
            self.filters_comments = [self.meta["filters"][k] for k in self.filters]

    def __repr__(self):
        name = self.filter_system
        if self.filters is None:
            name = self.filter_system + " at SVO"
        s = "Filter System: " + name

        return s


class Filter:

    def __init__(self, filter_name):

        self.filter_name = filter_name
        if self.filter_name in FILTER_DEFAULTS:
            self.filter_name = FILTER_DEFAULTS[filter_name]

        if "/" not in self.filter_name:
            raise ValueError("not a valid filter %s" % self.filter_name)

        self.filter_system, self.filter = self.filter_name.split("/")
        self.data_type = None
        self.wave_unit = None
        self.meta = None
        self.file_extension = None
        self.filename = None
        self.instrument = None

        self._update_atributes()
        self.path = self.get_path()

    def get_path(self):
        """ This just obtains the path to the filter  """
        relpath = urljoin("filter_systems", self.filter_system, self.filename)
        database = SpecDatabase()
        if self.filter_system in database.filter_systems:
            try:
                path = database.abspath(relpath)
            except (URLError, HTTPError):
                raise ValueError("filter %s is unknown" % self.filter_name)
        #elif self.filter_name in get_filter_names():
        else:
            download_file('http://svo2.cab.inta-csic.es/'
                          'theory/fps3/fps.php?ID={}'.format(self.filter_name),
                          os.path.join(database.rootdir, relpath))
            path = database.abspath(relpath)
       # else:
       #     raise ValueError("filter %s is unknown" % self.filter_name)

        return path

    def _update_atributes(self):
        database = SpecDatabase()
        if self.filter_system in database.filter_systems:
            system = FilterSystem(self.filter_system)
            self.wave_unit = system.wave_unit
            self.data_type = system.data_type

            self.file_extension = system.file_extension
            self.wave_unit = system.wave_unit
            self.instrument = system.instrument
            self.filename = self.filter + self.file_extension
            self.meta = system.meta

            del self.meta["filters"]

        else:
            self.filename = self.filter

    def __repr__(self):
        s = "Filter: " + self.filter_name
        return s


class ExtCurvesLibrary:
    def __init__(self, extcurvname):
        self.extcurvname = extcurvname
        database = SpecDatabase()
        if self.extcurvname not in database.extinction_curves:
            raise ValueError("Extinction curves not known")
        relpath = urljoin("extinction_curves", self.extcurvname, "index.yml")
        self.meta = database.get_yaml_contents(relpath)

        self.name = self.meta["name"]
        self.title = self.meta["title"]
        self.summary = self.meta["summary"]
        self.source = self.meta["source"]
        self.wave_unit = self.meta["wave_unit"]
        self.extinction_unit = self.meta["extinction_unit"]
        self.data_type = self.meta["data_type"]
        self.file_extension = self.meta["file_extension"]
        self.curves = list(self.meta["curves"].keys())
        self.curves_comments = [self.meta["curves"][k] for k in self.curves]

    def __repr__(self):
        s1 = "Extincion curves: " + self.name
        s2 = "Curves: " + str(self.curves)
        return '%s \n %s' % (s1, s2)


class ExtinctionCurve:

    def __init__(self, curve_name):
        self.curve_name = curve_name
        self.family, self.ext_curve = self.curve_name.split("/")
        ext_family = ExtCurvesLibrary(self.family)
        self.meta = ext_family.meta
        self.file_extension = self.meta["file_extension"]
        self.filename = self.ext_curve + self.file_extension
        self.path = self.get_path()

    def get_path(self):

        relpath = urljoin("extinction_curves", self.family, self.filename)
        database = SpecDatabase()
        path = database.abspath(relpath)

        return path

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
        flat_list = [f for f in filter_list if s in f]

    return flat_list
