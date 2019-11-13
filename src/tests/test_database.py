import pytest
import spyctra
import urllib
from astropy.table import Table


def test_database_location():
    """
    test if a valid url
    """
    url = spyctra.database.database_url()
    result = urllib.parse.urlparse(url)
    assert all([result.scheme, result.netloc, result.path])


database = spyctra.database.SpecDatabase()


class TestDatabase:

    def test_table(self):
        t = database.as_table
        assert isinstance(t, Table)


