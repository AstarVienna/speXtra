"""
Tests for class Spextrum
"""
import pytest
import inspect
import os

import numpy as np
import astropy.units as u
from astropy.constants import c
from synphot import SpectralElement, SourceSpectrum, units

import spextra
from spextra import Spextrum, Passband


# utility function
def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks"
    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))


MOCK_DIR = mock_dir()


def test_make_passband_no_file():
    with pytest.raises(FileNotFoundError) as e_info:
        pb = Passband.from_file(filename="blablabla")
        print(e_info)



class TestPassbandInstances:

    def test_alias(self):
        passband = Passband("W1")
        assert isinstance(passband, SpectralElement)

    def test_svo(self):
        passband = Passband("Paranal/HAWKI.Ks")
        assert isinstance(passband, SpectralElement)

    def test_database(self):
        passband = Passband("elt/micado/Y")
        assert isinstance(passband, SpectralElement)

    def test_filename(self):
        filter_file = os.path.join(MOCK_DIR, 'Y.dat')

        passband = Passband.from_file(filename=filter_file)
        assert isinstance(passband, SpectralElement)


class TestSpextrumInstances:
    """
    This class tests whether each method return the correct instances
    it also tests whether the methods are functional
    but it doesn't test for correct outputs see cls:TestSpextrum for that
    """

    sp = Spextrum("kc96/s0")

    def test_load(self, sp=sp):
        assert isinstance(sp, Spextrum)

    def test_sub_class(self):
        assert issubclass(Spextrum, SourceSpectrum)

    def test_redshift(self, sp=sp):

        sp2 = sp.redshift(z=1)
        assert isinstance(sp2, Spextrum)

    def test_add_emi_lines(self, sp=sp):
        sp2 = sp.add_emi_lines([5000, 6000], [10, 20], [1e-15, 2e-15])
        assert isinstance(sp2, Spextrum)

    def test_add_abs_lines(self, sp=sp):
        sp2 = sp.add_abs_lines([5000, 6000], [15, 20], [10, 12])
        assert isinstance(sp2, Spextrum)

    @pytest.mark.parametrize("unit", [u.mag, u.STmag, u.ABmag])
    def test_flat_spectrum(self, unit):
        sp = Spextrum.flat_spectrum(amplitude=10*unit)
        assert isinstance(sp, Spextrum)

    def test_mul_with_scalar(self, sp=sp):
        sp2 = sp * 2
        assert isinstance(sp2, Spextrum)

    def test_sum_spectra(self):
        sp1 = Spextrum("kc96/s0")
        sp2 = Spextrum("pickles/a0v")
        sp = sp1 + sp2
        assert isinstance(sp, Spextrum)

    def test_scale_to_magnitude(self, sp=sp):
        sp2 = sp.scale_to_magnitude(amplitude=15*u.ABmag, filter_curve="g")
        assert isinstance(sp2, Spextrum)

    def test_rebin_spectra(self, sp=sp):
        new_waves = np.linspace(np.min(sp.waveset),
                                np.max(sp.waveset),
                                100)
        sp2 = sp.rebin_spectra(new_waves=new_waves)
        assert isinstance(sp2, Spextrum)

    def test_get_magnitude(self):
        sp = Spextrum("pickles/a0v")
        mag = sp.get_magnitude("elt/micado/Y", system_name="AB")
        assert isinstance(mag, u.Quantity)

    def test_black_body_spectrum(self):
        sp = Spextrum.black_body_spectrum(filter_curve="g")
        assert isinstance(sp, Spextrum)

    def test_photons_in_range(self):
        sp = Spextrum.black_body_spectrum(filter_curve="g")
        counts = sp.photons_in_range(wmin=4000, wmax=5000)
        assert isinstance(counts, u.Quantity)

    def test_smooth(self):
        sp = Spextrum.black_body_spectrum(filter_curve="g")
        sp2 = sp.smooth(10*(u.m / u.s))
        assert isinstance(sp2, Spextrum)

    def test_redden(self):
        sp = Spextrum.black_body_spectrum(filter_curve="r")
        sp2 = sp.redden("calzetti/starburst", Ebv=0.1)
        assert isinstance(sp2, Spextrum)

    def test_deredden(self):
        sp = Spextrum.black_body_spectrum(filter_curve="F110W")
        sp2 = sp.redden("gordon/smc_bar", Ebv=0.1)
        assert isinstance(sp2, Spextrum)

    def testing_nesting(self):

        sp = Spextrum("kc96/s0").redshift(z=1).\
            scale_to_magnitude(amplitude=15*u.ABmag, filter_curve="g").\
            redden("calzetti/starburst", Ebv=0.1)
        assert isinstance(sp, Spextrum)





class TestSpextrum:

    def test_wrong_load(self):
        with pytest.raises(ValueError) as e_info:
            sp = Spextrum("kc96/wrong_name")

    @pytest.mark.parametrize("unit", [u.mag, u.STmag, u.ABmag])
    def test_ref_spectrum_is_right(self, unit):
        sp1 = Spextrum.flat_spectrum(amplitude=10*unit)
        sp2 = Spextrum.flat_spectrum(amplitude=11*unit)
        if unit == u.mag:
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

    @pytest.mark.parametrize("unit", [u.mag, u.ABmag, u.STmag])
    def test_correct_scaling(self, unit):
        sp1 = Spextrum("kc96/s0").scale_to_magnitude(amplitude=14*unit, filter_curve="r")
        sp2 = Spextrum("kc96/s0").scale_to_magnitude(amplitude=15*unit, filter_curve="r")

        flux1 = sp1(sp1.waveset[(sp1.waveset.value > 6231 - 200) &
                                (sp1.waveset.value < 6231 + 200)]).value
        flux2 = sp2(sp2.waveset[(sp2.waveset.value > 6231 - 200) &
                                (sp2.waveset.value < 6231 + 200)]).value

        mean = np.mean(flux1 / flux2)
        assert np.isclose(mean, 10**0.4)

    @pytest.mark.parametrize("filt", ["U", "B", "V", "R", "I", "J", "H", "Ks"])
    def test_vega2ab(self, filt):
        """
        test if the convertion between AB and Vega is correct
        conversions taken from:  http://www.astronomy.ohio-state.edu/~martini/usefuldata.html

        absolute tolerance set to 0.15 mag to account for filter differences

        Parameters
        ----------
        filt: str
            name of the filter
        """
        ab2vega = {"U":  0.79,   # magAB - magVega taken from
                   "B": -0.09,   #
                   "V":  0.02,
                   "R":  0.21,
                   "I":  0.45,
                   "J":  0.91,
                   "H":  1.39,
                   "Ks": 1.85}

        sp = Spextrum.flat_spectrum(amplitude=0*u.ABmag)

        magAB = sp.get_magnitude(filt, system_name="AB")
        magVega = sp.get_magnitude(filt, system_name="vega")

        diff = (magAB.value - magVega.value)

        assert np.isclose(diff, ab2vega[filt], atol=0.15)

    @pytest.mark.parametrize("spec_name", ["vega", "sun", "brown/NGC5992"])
    @pytest.mark.parametrize("filters", [["2MASS/2MASS.Ks", "elt/micado/Ks"],
                                         ["J", "elt/micado/J"],
                                         ["Generic/Bessell.B", "B"],
                                         ])
    def test_magnitudes(self, spec_name, filters):
        sp = Spextrum(spec_name)
        filter1 = filters[0]
        filter2 = filters[1]
        mag1 = sp.get_magnitude(filter_curve=filter1)
        mag2 = sp.get_magnitude(filter_curve=filter2)

        assert np.isclose(mag1.value, mag2.value, atol=0.3)

    def test_hz2angstrom(self):

        waves = np.array([1, 2, 3]) * u.Hz
        flux = np.array([1, 1, 2]) * units.FLAM

        sp = Spextrum.from_arrays(waves, flux)

        inwaves  = c.value / waves.value
        outwaves = sp.waveset.to(u.m).value

        assert np.isclose(inwaves[::-1], outwaves).all()





