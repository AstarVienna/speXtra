# -*- coding: utf-8 -*-
"""Pytest setup and global fixtures."""

from pathlib import Path

import pytest


@pytest.fixture(scope="package")
def mock_dir():
    return Path(__file__).parent / "mocks"
