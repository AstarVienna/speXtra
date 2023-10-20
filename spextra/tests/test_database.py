# -*- coding: utf-8 -*-
"""Tests for spextra database."""

import urllib
import pytest
import yaml

from spextra.database import Database
from spextra.libraries import SpecLibrary
from spextra.configuration import config


@pytest.fixture(scope="class", name="datacont")
def simple_yamldict(mock_dir):
    path = mock_dir / "index.yml"
    with path.open(encoding="utf-8") as file:
        return yaml.safe_load(file)


# @pytest.fixture(scope="module")
def mock_database():
    return Database()


# This doesn't work :(
# @pytest.fixture(scope="class")
# @pytest.mark.usefixtures("mock_database")
# def libraries(mock_database):
#     for library_name in mock_database.libraries:
#         yield library_name, SpecLibrary(library_name)


@pytest.fixture(scope="module")
def url_result():
    url = config.database_url
    return urllib.parse.urlparse(url)


@pytest.mark.parametrize("attribute", ["scheme", "netloc", "path"])
def test_database_location(url_result, attribute):
    """test if a valid url"""
    assert getattr(url_result, attribute)


class TestDatabase:
    @pytest.mark.parametrize("key", ["libraries",
                                     "extinction_curves",
                                     "filter_systems"])
    def test_keys(self, datacont, key):
        assert key in datacont
        assert len(datacont[key])


class TestSpecLibrary:
    @pytest.mark.parametrize("library_name", mock_database()["libraries"])
    def test_name(self, library_name):
        lib = SpecLibrary(library_name)
        assert lib.name == library_name

    def test_templates(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert "bulge" in lib.keys()

    @pytest.mark.parametrize("library_name", mock_database()["libraries"])
    @pytest.mark.parametrize("attribute", ["file_extension", "data_type"])
    def test_attr(self, library_name, attribute):
        lib = SpecLibrary(library_name)
        assert hasattr(lib, attribute)
