"""
Tests for class Spectrum
"""
import pytest
from spyctra import Spectrum, make_passband
from synphot import SpectralElement
import astropy.units as u

class TestPassbandInstances:
    def test_alias(self):
        passband = make_passband("W1")
        assert isinstance(passband, SpectralElement)

    def test_svo(self):
        passband = make_passband("Paranal/HAWKI.Ks")
        assert isinstance(passband, SpectralElement)

    def test_database(self):
        passband = make_passband("micado/Y")
        assert isinstance(passband, SpectralElement)

    def test_filename(self):
        filename = "/home/mverdugo/MyCodes/Spyctra/database/filter_systems/micado/Y.dat"
        passband = make_passband(filter_file=filename)
        assert isinstance(passband, SpectralElement)




class TestSpectrumInstances:

    def test_load(self):
        sp = Spectrum("kc96/s0")
        assert isinstance(sp, Spectrum)

    def test_redshift(self):
        sp = Spectrum("kc96/s0")
        sp2 = sp.redshift(z=1)
        assert isinstance(sp2, Spectrum)

    def test_add_emi_lines(self):
        sp = Spectrum("kc96/s0")
        sp2 = sp.add_emi_lines(5000, 1e-15, 10 )
        assert isinstance(sp2, Spectrum)

    def test_add_abs_lines(self):
        sp = Spectrum("kc96/s0")
        sp2 = sp.add_abs_lines(5000, 15, 10 )
        assert isinstance(sp2, Spectrum)

    def test_ref_spectrum(self):
        sp = Spectrum.ref_spectrum(mag=10, system_name="vega")
        assert isinstance(sp, Spectrum)

    def test_scale_to_magnitude(self):
        sp = Spectrum("kc96/s0")
        sp2 = sp.scale_to_magnitude(amplitude=13*u.AB, filter_name="g")
        assert isinstance(sp2, Spectrum)

class TestSpectrum:

    def test_wrong_load(self):
        with pytest.raises(ValueError) as e_info:
            sp = Spectrum("kc96/wrong_name")

    def test_units(self):
        pass

    def test_scale_to_magnitude(self):
        pass





