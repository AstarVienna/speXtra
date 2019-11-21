"""
Tests for class Spectrum
"""
import pytest
from spyctra import Spectrum, make_passband
from synphot import SpectralElement, SourceSpectrum, units
import astropy.units as u
import numpy as np


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
    """
    This class tests whether each method return the correct instances
    it also tests whether the methods are functional
    but it doesn't test for correct outputs see cls:TestSpectrum for that
    """

    sp = Spectrum("kc96/s0")

    def test_load(self, sp=sp):
        assert isinstance(sp, Spectrum)

    def test_sub_class(self):
        assert issubclass(Spectrum, SourceSpectrum)

    def test_redshift(self, sp=sp):
        sp2 = sp.redshift(z=1)
        assert isinstance(sp2, Spectrum)

    def test_add_emi_lines(self, sp=sp):
        sp2 = sp.add_emi_lines(5000, 1e-15, 10 )
        assert isinstance(sp2, Spectrum)

    def test_add_abs_lines(self, sp=sp):
        sp2 = sp.add_abs_lines(5000, 15, 10 )
        assert isinstance(sp2, Spectrum)

    @pytest.mark.parametrize("system_name", ["ab", "st", "vega"])
    def test_ref_spectrum(self, system_name):
        sp = Spectrum.ref_spectrum(mag=10, system_name=system_name)
        assert isinstance(sp, Spectrum)

    def test_mul_with_scalar(self, sp=sp):
        sp2 = sp * 2
        assert isinstance(sp2, Spectrum)

    def test_sum_spectra(self):
        sp1 = Spectrum("kc96/s0")
        sp2 = Spectrum("pickles/a0v")
        sp = sp1 + sp2
        assert isinstance(sp, Spectrum)

    def test_scale_to_magnitude(self, sp=sp):
        sp2 = sp.scale_to_magnitude(amplitude=13*u.ABmag, filter_name="g")
        assert isinstance(sp2, Spectrum)

    def test_rebin_spectra(self, sp=sp):
        new_waves = np.linspace(np.min(sp.waveset),
                                np.max(sp.waveset),
                                100)
        sp2 = sp.rebin_spectra(new_waves=new_waves)
        assert isinstance(sp2, Spectrum)

    def test_get_magnitude(self, sp=sp):
        mag = sp.get_magnitude("g", system_name="AB")
        assert isinstance(mag, u.Quantity)


class TestSpectrum:

    def test_wrong_load(self):
        with pytest.raises(ValueError) as e_info:
            sp = Spectrum("kc96/wrong_name")

    @pytest.mark.parametrize("system_name", ["ab", "st", "vega"])
    def test_ref_spectrum_is_right(self, system_name):
        sp1 = Spectrum.ref_spectrum(mag=10, system_name=system_name)
        sp2 = Spectrum.ref_spectrum(mag=11, system_name=system_name)
        if system_name == "vega":
            flux1 = sp1(sp1.waveset[(sp1.waveset.value > 7000 - 200) &
                                    (sp1.waveset.value < 7000 + 200)]).value
            flux2 = sp2(sp2.waveset[(sp2.waveset.value > 7000 - 200) &
                                    (sp2.waveset.value < 7000 + 200)]).value

        else:
            waves = np.arange(1000, 1e4, 1) * u.AA
            flux1 = sp1(waves)
            flux2 = sp2(waves)

        mean = np.mean(flux1 / flux2)
        assert np.isclose(mean, 10**0.4)

    @pytest.mark.parametrize("units", [u.mag, u.ABmag, u.STmag])
    def test_scaling_is_right(self, units):
        sp1 = Spectrum("kc96/s0").scale_to_magnitude(amplitude=14*units, filter_name="r")
        sp2 = Spectrum("kc96/s0").scale_to_magnitude(amplitude=15*units, filter_name="r")

        flux1 = sp1(sp1.waveset[(sp1.waveset.value > 6231 - 200) &
                                (sp1.waveset.value < 6231 + 200)]).value
        flux2 = sp2(sp2.waveset[(sp2.waveset.value > 6231 - 200) &
                                (sp2.waveset.value < 6231 + 200)]).value

        mean = np.mean(flux1 / flux2)
        assert np.isclose(mean, 10**0.4)

    def test_units(self):
        pass

    def test_scale_to_magnitude(self):
        pass





