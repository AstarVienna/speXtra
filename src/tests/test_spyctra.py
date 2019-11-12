"""
Tests for class Spectrum
"""
import pytest
from spyctra import Spectrum


class TestSpectrum:

    def test_load(self):
        sp = Spectrum.load("kc96/s0")
        assert isinstance(sp, Spectrum)

    def test_wrong_load(self):
        with pytest.raises(ValueError) as e_info:
            sp = Spectrum.load("kc96/wrong_name")

    def test_redshift(self):
        sp = Spectrum.load("kc96/s0")
        sp2 = sp.redshift(z=1)
        assert isinstance(sp2, Spectrum)