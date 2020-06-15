# -*- coding: utf-8 -*-
"""
Database for speXtra
"""
import inspect
from posixpath import join as urljoin
import os
import yaml
from urllib.error import URLError, HTTPError

import tynt

from .utils import get_rootdir, database_url, download_file

# Configurations

__all__ = ["SpecDatabase", "SpecLibrary", "SpectralTemplate",
           "ExtinctionCurve", "Filter"]

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

        self._checked_rootdir = None
        self.rootdir = rootdir
        self.remote_root = remote_root

        self.url = remote_root + "index.yml"
        self.contents = self.get_yaml_contents("index.yml")
        self.libraries = [lib for lib in self.contents["libraries"]]
        self.extinction_curves = [ext for ext in self.contents["extinction_curves"]]
        self.filter_systems = [filt for filt in self.contents["filter_systems"]]

    def abspath(self, relpath):
        """Return absolute path to file or directory, ensuring that it exists.
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
        path

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
        libs = "libraries:" + ' ' + str(self.library_names)
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
        self.location = urljoin(database_url(), "libraries",
                                self.name, "index.yml")
        self.data = self.get_data()
        self.library_name = self.data["library_name"]
        self.title = self.data["title"]
        self.type = self.data["type"]
        self.summary = self.data["summary"]
        self.reference = self.data["reference"]
        self.link = self.data["link"]
        self.spectral_coverage = self.data["spectral_coverage"]
        self.resolution = self.data["resolution"]
        self.wave_unit = self.data["wave_unit"]
        self.flux_unit = self.data["flux_unit"]
        self.wave_column_name = self.data["wave_column_name"]
        self.flux_column_name = self.data["flux_column_name"]
        self.data_type = self.data["data_type"]
        self.file_extension = self.data["file_extension"]
        self.templates = list(self.data["templates"].keys())
        self.template_comments = [self.data["templates"][k] for k in self.templates]

    def get_data(self):
        relpath = urljoin("libraries", self.name, "index.yml")
        database = SpecDatabase(get_rootdir(), database_url())

        return database.get_yaml_contents(relpath)

    def dump(self):
        """
        Nicely dump the contents of the library

        Returns
        -------

        """
        print(yaml.dump(self.data,
                        indent=4, sort_keys=False, default_flow_style=False))

    def __repr__(self):
        description = "Spectral Library: " + self.name + " " + self.title
        spec_cov = "spectral coverage: " + str(self.spectral_coverage)
        units = "wave_unit: " + self.wave_unit + "  flux_unit: " + self.flux_unit
        templates = "Templates: " + str(self.templates)

        return ' %s \n %s \n %s \n %s' % (description, spec_cov, units, templates)


class SpectralTemplate:

    def __init__(self, template):
        self.library_name, self.template_name = template.split("/")
        library = SpecLibrary(self.library_name)
        self.resolution = library.resolution
        self.wave_unit = library.wave_unit
        self.flux_unit = library.flux_unit
        self.wave_column_name = library.wave_column_name
        self.flux_column_name = library.flux_column_name
        self.data_type = library.data_type
        self.file_extension = library.file_extension
        self.filename = self.template_name + self.file_extension
        self.path = self.get_data()

    def get_path(self):
        relpath = urljoin("libraries", self.library_name, self.filename)
        database = SpecDatabase(get_rootdir(), database_url())

        return database.abspath(relpath)

    def __repr__(self):
        s1 = "Spectral template: " + self.template
        return s1


class FilterSystem:

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
        self._update_atributes()

    def get_data(self):
        data = None
        relpath = urljoin("filter_systems", self.filter_system, "index.yml")
        database = SpecDatabase(get_rootdir(), database_url())
        if self.filter_system in database.filter_systems:
            data = database.get_yaml_contents(relpath)

        return data

    def _update_atributes(self):
        data = self.get_data()
        if data is not None:
            self.filter_system = data["filter_system"]
            self.instrument = data["instrument"]
            self.title = data["title"]
            self.author = data["author"]
            self.source = data["source"]
            self.spectral_coverage = data["spectral_coverage"]
            self.wave_unit = data["wave_unit"]
            self.file_extension = data["file_extension"]
            self.filters = list(data["filters"].keys())
            self.filters_comments = [data["filters"][k] for k in self.filters]

    def __repr__(self):
        name = self.filter_system
        if self.filters is None:
            name = self.filter_system + " at SVO"
        s = "Filter System: " + name

        return s

#




class Filter:

    def __init__(self, filter_name):

        self.filter_name = filter_name
        try:
            self.filter_system, self.filter = self.filter_name.split("/")

        except ValueError:
            if filter_name in FILTER_DEFAULTS:
                self.filter_name = FILTER_DEFAULTS[filter_name]
                self.filter_system, self.filter = self.filter_name.split("/")

        self.data = self.get_data()

    def get_data(self):
        relpath = urljoin("filters", self.filter_system, self.filter)
        try:
            database = SpecDatabase(get_rootdir(), database_url())
            path = database.abspath(relpath)
        except (URLError, HTTPError):
            path = download_file('http://svo2.cab.inta-csic.es/'
                                 'theory/fps3/fps.php?ID={}'.format(self.filter_name),
                                 os.path.join(get_rootdir(), relpath))

        return path

    def __repr__(self):
        s = "Filter: " + self.filter_name
        return s


class ExtCurvesLibraries:
    def __init__(self, extcurvname):
        self.extcurvname = extcurvname
        pass

    def get_data(self):
        relpath = urljoin("extinction_curves", self.extcurvname)
        database = SpecDatabase(get_rootdir(), database_url())
        data = database.get_yaml_contents(relpath)
        return data


class ExtinctionCurve:

    def __init__(self, curve_name):
        self.curve_name = curve_name
        self.family, self.ext_curve = self.curve_name.split("/")
        self.path = self.getdata()
        self.meta = None

    def get_path(self):
        relpath = urljoin("extinction_curves", self.family, self.ext_curve)
        database = SpecDatabase(get_rootdir(), database_url())
        path = database.get_yaml_contents(relpath)

        return path

    def __repr__(self):
        s = "Extinction curve: " + self.curve_name
        return s








# This is based on scopesim.effects.ter_curves_utils.py


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
