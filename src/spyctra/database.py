# -*- coding: utf-8 -*-
"""
Database for spyctra
"""
import shutil
from posixpath import join as urljoin
import urllib

import yaml
from astropy.utils.data import download_file
from astropy.table import Table


def is_url(url):
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
        self.library_names = [l for l in self.contents["library_names"]]

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

        path = urljoin("libraries", library_name, "contents.yml")
        return self._get_contents(path)

    def display_library(self, library_name):
        """
        try to nicely display the contents of library

        """
        print(yaml.dump(self.get_library(library_name),
                        indent=4, sort_keys=False, default_flow_style=False))

    @property
    def as_table(self):
        """
        make a summary of the database properties

        Returns
        -------
        an astropy.table with the database contents
        """

        column_names = ["library_name", "title", "type", "resolution", "spectral_coverage", "templates"]
        library_names = self.library_names
        titles = []
        types = []
        resolution = []
        spectral_coverage = []
        templates = []

        for lib in self.library_names:
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
        Try to nicely display the whole database
        Returns
        -------

        """
        print(yaml.dump(self.as_dict,
                        indent=4, sort_keys=False, default_flow_style=False))

    def browse(self, keys):
        """

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









