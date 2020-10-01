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
    """
    def __init__(self, filename):
        self.filename = filename
        with open(filename) as f:
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
        self.filename = self.abspath("index.yml")

        super().__init__(filename=self.filename)

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
        self.library_name = library_name
        self.relpath = os.path.join("libraries", library_name, "index.yml")
        self.path = database.abspath(self.relpath)
        super().__init__(filename=self.path)

    def __repr__(self):
        description = "Spectral Library: " + self.library_name + " " + self.title
        spec_cov = "spectral coverage: " + str(self.spectral_coverage)
        units = "wave_unit: " + self.wave_unit + "  flux_unit: " + self.flux_unit
        templates = "Templates: " + str(self.templates)

        return ' %s \n %s \n %s \n %s' % (description, spec_cov, units, templates)


class FilterSystem(DataContainer):

    def __init__(self, filter_system):
        filter_system_path = os.path.join("filter_systems", filter_system)
        database = Database()
        if filter_system_path not in database.pathlist :
            raise ValueError("filter system '%s' is not in the database" % filter_system)

        self.relpath = os.path.join(filter_system_path, "index.yml")
        self.path = database.abspath(self.relpath)
        super().__init__(filename=self.path)


class ExtCurvesLibrary(DataContainer):

    def __init__(self, extinction_curve):
        database = Database()
        if extinction_curve not in database.extinction_curves:
            raise ValueError("extinction curve '%s' is not in the database" % extinction_curve)

        self.relpath = os.path.join("filter_systems", extinction_curve, "index.yml")
        self.path = database.abspath(self.relpath)
        super().__init__(filename=self.path)


class SpectrumContainer(SpecLibrary):
    """
    A holder of spectral template information, including template characteristics and location
    """

    def __init__(self, template):
        library = template.split("/")[:-1]
        self.template_name = template.split("/").pop()
        if self.template_name not in library.templates.keys():
            raise ValueError("Template not in library", self.template_name)

        super.__init__(filename=library)

        self.filename = self.template_name + self.file_extension
        self.path = self.get_path()

    def get_path(self):
        database = Database()
        relpath = urljoin(self.relpath, self.filename)
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


class FilterContainer:

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
        database = Database()
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
        database = Database()
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




class ExtCurveContainer:

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
        database = Database()
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
        flat_list = [f for f in filter_list if system in f]

    return flat_list
