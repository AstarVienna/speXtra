# -*- coding: utf-8 -*-
"""Tests for spextra database."""

import urllib
import pytest

from spextra.database import DataContainer, Database, SpecLibrary
from spextra.configuration import config


@pytest.fixture(scope="class", name="datacont")
@pytest.mark.usefixtures("mock_dir")
def simple_datacontainer(mock_dir):
    return DataContainer(filename=mock_dir / "index.yml")


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
    url = config.get_database_url()
    return urllib.parse.urlparse(url)


@pytest.mark.usefixtures("url_result")
@pytest.mark.parametrize("attribute", ["scheme", "netloc", "path"])
def test_database_location(url_result, attribute):
    """test if a valid url"""
    assert getattr(url_result, attribute)


class TestDataContainer:
    @pytest.mark.usefixtures("datacont")
    @pytest.mark.parametrize("attribute", ["libraries",
                                           "extinction_curves",
                                           "filter_systems"])
    def test_attr_creation(self, datacont, attribute):
        assert hasattr(datacont, attribute)


class TestSpecLibrary:
    @pytest.mark.parametrize("library_name", mock_database().libraries)
    def test_name(self, library_name):
        lib = SpecLibrary(library_name)
        assert lib.name == library_name

    def test_templates(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert "bulge" in lib.items.keys()

    @pytest.mark.parametrize("library_name", mock_database().libraries)
    @pytest.mark.parametrize("attribute", ["file_extension", "data_type"])
    def test_attr(self, library_name, attribute):
        lib = SpecLibrary(library_name)
        assert hasattr(lib, attribute)
