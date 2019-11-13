"""
Tests for class Spectrum
"""
import pytest
from spyctra import Spectrum


class TestSpectrumInstances:

    def test_load(self):
        sp = Spectrum("kc96/s0")
        assert isinstance(sp, Spectrum)

    def test_wrong_load(self):
        with pytest.raises(ValueError) as e_info:
            sp = Spectrum("kc96/wrong_name")

    def test_redshift(self):
        sp = Spectrum("kc96/s0")
        sp2 = sp.redshift(z=1)
        assert isinstance(sp2, Spectrum)

    def test_zero_mag_spectrum(self):
        sp = Spectrum.zero_mag_spectrum(system_name="AB")
        assert isinstance(sp, Spectrum)



