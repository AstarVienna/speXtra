# -*- coding: utf-8 -*-
"""
Database for speXtra
"""
import shutil
from posixpath import join as urljoin
import urllib
import os
import inspect

import yaml
from astropy.utils.data import download_file
from astropy.table import Table

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

# default filter definitions now in data/default_filters.yml,
# including now also GALEX, Spitzer, WISE HST, Y-band and GAIA G-band
with open(os.path.join(__data_dir__,  "default_filters.yml")) as filter_file:
    FILTER_DEFAULTS = yaml.safe_load(filter_file)


def is_url(url):  # TODO: Move to .utils
    """
    Checks that a given URL is reachable.
    Depending of the configuration of the server it might return True even if the page doesn't exist.

    Parameters
    -----------
    url: A URL

    Returns
    -------
    Boolean
    """
    try:
        request = urllib.request.Request(url)
        request.get_method = lambda: 'HEAD'
        urllib.request.urlopen(request)
        output = True
    except urllib.error.URLError:
        output = False
    except ValueError:
        output = False
    finally:
        return output


def database_url():
    """
    TODO: it should read it from a file
    Returns
    -------
    the database_location
    """
    loc = "https://homepage.univie.ac.at/miguel.verdugo/database/"
    try:
        assert is_url(loc)
    except AssertionError as error:
        print(error)
        print("Database address not reachable", loc)
    return loc


class SpecDatabase:

    def __init__(self):
        self.url = database_url()
        self.library_names = [lib for lib in self.contents["library_names"]]
        self.extinction_curves = [ext for ext in self.contents["extinction_curves"]]
        self.filter_systems = [filt for filt in self.contents["filter_systems"]]

    @property
    def contents(self):
        return self._get_contents("index.yml")

    def get_library(self, library_name):
        """
        get the contents of a particular library
        Parameters
        ----------
        library_name

        Returns
        -------
        a dictionary with the library contents
        """

        if library_name not in self.library_names:
            raise ValueError(library_name, "library not found")

        path = urljoin("libraries", library_name, "index.yml")
        return self._get_contents(path)

    def get_extinction_curves(self, extinction_name):

        if extinction_name not in self.extinction_curves:
            raise ValueError(extinction_name, "extinction curves not found")
        pass

    def get_filter_system(self, filter_system):

        if filter_system not in self.filter_systems:
            raise ValueError(filter_system, "filter system not found")

        path = urljoin("filter_systems", filter_system, "index.yml")
        return self._get_contents(path)

    def print_library(self, library_name):
        """
        just try to nicely print the contents of a library
        """
        print(yaml.dump(self.get_library(library_name),
                        indent=4, sort_keys=False, default_flow_style=False))

    def print_filter_system(self, filter_system):
        print(yaml.dump(self.get_filter_system(filter_system),
                        indent=4, sort_keys=False, default_flow_style=False))

    @property
    def libraries_as_table(self):
        """
        make a summary of the libraries properties

        Returns
        -------
        an astropy.table.Table with the database contents
        """

        column_names = ["library_name", "title", "type", "resolution", "spectral_coverage", "templates"]
        library_names = self.library_names
        titles = []
        types = []
        resolution = []
        spectral_coverage = []
        templates = []

        for lib in library_names:
            contents = self.get_library(lib)
            titles.append(contents["title"])
            types.append(contents["type"])
            resolution.append(contents["resolution"])
            spectral_coverage.append(contents["spectral_coverage"])
            templates.append([t for t in contents["templates"]])

        data = [library_names, titles, types, resolution, spectral_coverage, templates]
        table = Table(names=column_names, data=data)

        return table

    @property
    def filters_as_table(self):
        """
        make a summary of the filters available properties
        Returns
        -------
        an astropy.table.Table
        """
        column_names = ["filter_system", "instrument", "title", "spectral_coverage", "filters"]
        filter_systems = self.filter_systems
        instruments = []
        titles = []
        spectral_coverage = []
        filters = []
        for filt in filter_systems:
            contents = self.get_filter_system(filt)
            instruments.append(contents["instrument"])
            titles.append(contents["title"])
            spectral_coverage.append(contents["spectral_coverage"])
            filters.append([f for f in contents["filters"]])

        data = [filter_systems, instruments, titles, spectral_coverage, filters]
        table = Table(names=column_names, data=data)

        return table

    @property
    def as_dict(self):
        """
        Represent the whole database as a dictionary

        Returns
        -------

        """
        database = {}
        for lib in self.library_names:
            contents = self.get_library(lib)
            database[lib] = contents

        return database

    def display(self):
        """
        just try to nicely print the contents of the whole database

        Returns
        -------

        """
        print(yaml.dump(self.as_dict,
                        indent=4, sort_keys=False, default_flow_style=False))

    def browse(self, keys):
        """
        TODO: Implement
        Parameters
        ----------
        keys

        Returns
        -------
        libraries and templates that fulfill the criteria (type, spectral coverage, resolution etc
        """
        pass

    def _get_contents(self, path=""):
        """
        read a yaml file from a relative url

        Parameters
        ----------
        path

        Returns
        -------
        dict with the contents of the yaml file
        """

        url = urljoin(self.url, path)
        try:
            assert is_url(url)
        except urllib.error.URLError:
            print(url, "library not reachable")
        else:
            filename = download_file(url, cache=True)
            with open(filename) as f:
                data = yaml.safe_load(f)

        return data


def get_template(template, path=None):
    """
    Parameters
    ----------
    template: the name of the template, specified as library_name/template_name
    path: the path were the downloaded template will be downloaded (optional)

    Returns
    -------
    a file
    a dictionary with the main template attributes

    """
    newfile = None
    database = SpecDatabase()
    library_name, template_name = template.split("/")
    lib_data = database.get_library(library_name)
    template_meta = {"resolution": lib_data["resolution"],
                     "wave_unit": lib_data["wave_unit"],
                     "flux_unit": lib_data["flux_unit"],
                     "wave_column_name": lib_data["wave_column_name"],
                     "flux_column_name": lib_data["flux_column_name"],
                     "data_type": lib_data["data_type"],
                     "file_extension": lib_data["file_extension"]}

    try:
        assert template_name in lib_data["templates"]
    except AssertionError as error:
        print(error)
        print(template_name, "not found")
    else:
        filename = template_name + template_meta["file_extension"]
        url = urljoin(database.url, "libraries/", library_name, filename)
        newfile = download_file(url, cache=True)
        if path is not None:
            file = shutil.copy2(newfile, path)
            print(filename)

    return newfile, template_meta


def get_filter(filter_name):
    """
    get filter from the database. It will try first to download from the speXtra database
    and then try to get it from the SVO service

    TODO: Probably needs a refactoring

    Parameters
    ----------
    filter_name: str with the following format "instrument/filter"

    Returns
    -------
    a path and metadata
    """
    database = SpecDatabase()

    try:
        filter_system, filt = filter_name.split("/")
        filter_data = database.get_filter_system(filter_system)
        filter_meta = {"wave_unit": filter_data["wave_unit"],
                       "file_extension": filter_data["file_extension"],
                       "data_type": filter_data["data_type"]}

        filename = filt + filter_meta["file_extension"]
        url = urljoin(database.url, "filter_systems/", filter_system, filename)
        path = download_file(url, cache=True)
    except ValueError:
        if filter_name in FILTER_DEFAULTS:
            filter_name = FILTER_DEFAULTS[filter_name]

        path = download_file('http://svo2.cab.inta-csic.es/'
                             'theory/fps3/fps.php?ID={}'.format(filter_name),
                             cache=True)
        filter_meta = None
        with open(path) as f:
            if "Filter not found" in f.read():
                path = None
                raise ValueError("Filter not found")

    return path, filter_meta


def get_extinction_curve(curve_name):
    """
    Download a extinction curve from the database

    Parameters
    ----------
    curve_name: extinction curve, e.g.

    Returns
    -------

    """
    pass


# This is based on scopesim.effects.ter_curves_utils.py

