"""
Tests for spextra database
"""
import urllib
import inspect
import os

import spextra
from spextra.database import SpecLibrary, SpecDatabase
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


database = spextra.database.SpecDatabase()


class TestSpecLibrary:

    def test_name(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert lib.name == name

    def test_location(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert is_url(lib.location)

    def test_data(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert isinstance(lib.meta, dict)

    def test_templates(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert "bulge" in lib.templates


class TestSpecDatabase:
    def test_remote_root(self):
        database = SpecDatabase(remote_root="http://www.google.com")
        assert database.remote_root == "http://www.google.com/"

    def test_rootdir(self):
        database = SpecDatabase(rootdir=MOCK_DIR)
        assert "mocks" in database.rootdir

    def test_abspath(self):
        database = SpecDatabase(rootdir=MOCK_DIR)
        assert os.path.exists("mocks/index.yml")

    def test_get_yaml_contents(self):
        database = SpecDatabase(rootdir=MOCK_DIR)
        data = database.get_yaml_contents("index.yml")
        assert isinstance(data, dict)
    
