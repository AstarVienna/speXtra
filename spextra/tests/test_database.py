import spextra
import urllib
from astropy.table import Table
from spextra.database import SpecLibrary, SpecDatabase, is_url

def test_database_location():
    """
    test if a valid url
    """
    url = spextra.database.database_url()
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

    def test_box(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert isinstance(lib.fields, dict)

    def test_fields(self):
        name = "kc96"
        lib = SpecLibrary(name)
        assert isinstance(lib.fields.templates, dict)








class TestDatabase:

    def test_table(self):
        t = database.as_table()
        assert isinstance(t, Table)
