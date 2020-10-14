"""
Tests for spextra database
"""
import urllib
import inspect
import os
import pytest

import spextra
from spextra.database import DataContainer, Database, Library, SpecLibrary, FilterSystem, ExtCurvesLibrary
from spextra.utils import is_url


def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks"
    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))


MOCK_DIR = mock_dir()


def test_database_location():
    """
    test if a valid url
    """
    url = spextra.utils.database_url()
    result = urllib.parse.urlparse(url)
    assert all([result.scheme, result.netloc, result.path])


class TestDataContainer:

    def test_attr_creation(self):
        dc = DataContainer(filename=os.path.join(MOCK_DIR, "index.yml"))

        assert hasattr(dc, "libraries")
        assert hasattr(dc, "extinction_curves")
        assert hasattr(dc, "filter_systems")


class TestDatabase:
    def test_remote_root(self):
        database = Database(remote_root="http://www.google.com")
        assert database.remote_root == "http://www.google.com/"

    def test_rootdir(self):
        database = Database(rootdir=MOCK_DIR)
        assert "mocks" in database.rootdir

    def test_meta(self):
        database = Database(rootdir=MOCK_DIR)
        data = database.meta
        assert isinstance(data, dict)


class TestLibrary:
    pass


class TestSpecLibrary:

    database = Database()
    libraries = database.libraries
    needed_attr = ["file_extension", "data_type"]   # list here the mandatory attributes TODO: Fill that

    @pytest.mark.parametrize("library_name", libraries)
    def test_name(self, library_name):

        lib = SpecLibrary(library_name)
        assert lib.library_name == library_name

    def test_templates(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert "bulge" in lib.template_names

    @pytest.mark.parametrize("library_name", libraries)
    @pytest.mark.parametrize("attribute", needed_attr)
    def test_attr(self, library_name, attribute):
        lib = SpecLibrary(library_name)
        assert hasattr(lib, attribute)


