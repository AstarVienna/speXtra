# -*- coding: utf-8 -*-
"""Pytest setup and global fixtures."""

import pytest

from . import PATH_HERE

@pytest.fixture(scope="package")
def mock_dir():
    return PATH_HERE / "mocks"
