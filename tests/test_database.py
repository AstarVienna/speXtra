import pytest
import spyctra
from urllib.request import urlopen

print(spyctra.database_location)

loc = 'http://www.google.cm'


def test_url():
    code = urlopen(loc).code
    assert code < 400

