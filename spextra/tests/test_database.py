# -*- coding: utf-8 -*-
"""Tests for spextra database."""

import urllib
import pytest
import yaml

from spextra.database import Database, DefaultData
from spextra.libraries import SpecLibrary
from spextra.configuration import config

from . import PATH_HERE

path = PATH_HERE / "mocks/index.yml"
with path.open(encoding="utf-8") as file:
    datacont = yaml.safe_load(file)


@pytest.fixture(scope="class")
def mock_database():
    return datacont


@pytest.fixture(scope="class")
def mock_defaultdata():
    return DefaultData.from_yamls()


@pytest.fixture(scope="module")
def url_result():
    url = config.database_url
    return urllib.parse.urlparse(url)


@pytest.mark.parametrize("attribute", ["scheme", "netloc", "path"])
def test_database_url(url_result, attribute):
    """test if a valid url"""
    assert getattr(url_result, attribute)


class TestDatabase:
    @pytest.mark.parametrize("key", ["libraries",
                                     "extinction_curves",
                                     "filter_systems"])
    def test_keys(self, key, mock_database):
        assert key in mock_database
        assert len(mock_database[key])


class TestDefaultData:
    @pytest.mark.parametrize("attr", ["filters", "extcurves", "spectra"])
    def test_attrs(self, attr, mock_defaultdata):
        assert hasattr(mock_defaultdata, attr)


class TestSpecLibrary:
    @pytest.mark.parametrize("library_name", datacont["libraries"])
    def test_name(self, library_name):
        lib = SpecLibrary(library_name)
        assert lib.name == library_name

    def test_templates(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert "bulge" in lib.keys()

    @pytest.mark.parametrize("library_name", datacont["libraries"])
    @pytest.mark.parametrize("attribute", ["file_extension", "data_type"])
    def test_attr(self, library_name, attribute):
        lib = SpecLibrary(library_name)
        assert hasattr(lib, attribute)
