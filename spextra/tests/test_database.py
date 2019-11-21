import spextra
import urllib
from astropy.table import Table


def test_database_location():
    """
    test if a valid url
    """
    url = spextra.database.database_url()
    result = urllib.parse.urlparse(url)
    assert all([result.scheme, result.netloc, result.path])


database = spextra.database.SpecDatabase()


class TestDatabase:

    def test_table(self):
        t = database.libraries_as_table
        assert isinstance(t, Table)
