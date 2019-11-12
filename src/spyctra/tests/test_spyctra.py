"""
Tests for class Spectrum
"""

import pytest
from spyctra import Spectrum


class TestSpectrum:

    def test_load(self):
        sp = Spectrum.load("kc98/s0")
        assert isinstance(Spectrum, sp)



