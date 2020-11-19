# -*- coding: utf-8 -*-
"""
speXtra: A python tool to manage and manipulate astronomical spectra
"""
import numbers
import warnings

import numpy as np
from scipy.ndimage import gaussian_filter1d

import astropy.units as u
import astropy.constants as constants
from astropy.modeling.models import Scale

from synphot import (units, SourceSpectrum, SpectralElement, Observation,
                     BaseUnitlessSpectrum, ReddeningLaw, utils)
from synphot.models import (Empirical1D, GaussianFlux1D, Gaussian1D, Box1D, ConstFlux1D, BlackBody1D, PowerLawFlux1D)
from synphot.specio import read_ascii_spec, read_fits_spec
from synphot import exceptions

from .database import SpectrumContainer, FilterContainer, ExtCurveContainer, DefaultData
from .utils import download_svo_filter


__all__ = ["Spextrum", "Passband", "ExtinctionCurve",  "get_vega_spectrum"]

speed_of_light = constants.c

DEFAULT_FILTERS = DefaultData().filters
DEFAULT_SPECTRA = DefaultData().spectra
DEFAULT_CURVES = DefaultData().extcurves


class Passband(SpectralElement, FilterContainer):
    """
    This should be the holder of all information and operations related to the filters
    including path, etc.

    TODO: Implement proper __add__, __sub__, __mul__, etc

    """

    def __init__(self, filter_name=None, modelclass=None, **kwargs):

        if filter_name is not None:
            if filter_name in DEFAULT_FILTERS.keys():
                filter_name = DEFAULT_FILTERS[filter_name]

            try:
                FilterContainer.__init__(self, filter_name)
                meta, wave, trans = self._loader()
                SpectralElement.__init__(self, Empirical1D, points=wave, lookup_table=trans, meta=meta)
            except ValueError as e1:
                try:  # try to download it from SVO
                    meta, wave, trans = self._from_svo(filter_name)
                    SpectralElement.__init__(self, Empirical1D, points=wave, lookup_table=trans,
                                             meta=meta)
                except ValueError as e2:
                    print("filter %s doesn't exist" % filter_name)

        elif modelclass is not None:
            SpectralElement.__init__(self, modelclass)
        else:
            raise ValueError("please define a filter")

    def _loader(self):
        """
        Load a filter from the database

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """

        try:  # it should also try to read it from the file directly
            self.wave_unit = units.validate_unit(self.wave_unit)
        except exceptions.SynphotError:
            self.wave_unit = u.AA

        # make try and except here to catch most problems
        if self.data_type == "fits":
            meta, lam, trans = read_fits_spec(self.filename, ext=1,
                                              wave_unit=self.wave_unit,
                                              wave_col=self.wave_column_name, flux_col=self.trans_column_name)
        elif self.data_type == "ascii":
            meta, lam, trans = read_ascii_spec(self.filename,
                                               wave_unit=self.wave_unit, flux_unit=self._internal_flux_unit,
                                             #  wave_col=self.wave_column_name, flux_col=self.trans_column_name,
                                               )

        return meta, lam, trans

    @classmethod
    def from_vectors(cls, waves, trans, meta=None, wave_unit=u.AA):
        """
        Create a ``Passband`` directly from from vectos (lists, numpy.arrays, etc)
        Parameters
        ----------
        waves: list-like
        trans: list-like
        meta: dictionary containing the metadata
        wave_unit: u.Quantity, defaulted to angstroms

        Returns
        -------
        Passband
        """
        if isinstance(waves, u.Quantity) is False:
            waves = waves*wave_unit

        modelclass = SpectralElement(Empirical1D, points=waves, lookup_table=trans, meta=meta)

        return cls(modelclass=modelclass)

    @classmethod
    def from_file(cls, filename, **kwargs):

        modelclass = SpectralElement.from_file(filename, **kwargs)

        return cls(modelclass=modelclass)

    @classmethod
    def gaussian(cls, center, fwhm, peak, **kwargs):
        """
        Creates a filter with a gaussian shape with given user parameters
        Returns
        -------

        """
        if isinstance(center, u.Quantity) is False:
            center = center*u.AA
        if isinstance(fwhm, u.Quantity) is False:
            fwhm = fwhm*u.AA

        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

        modelclass = SpectralElement(Gaussian1D, amplitude=peak, mean=center, stddev=sigma, **kwargs)

        return cls(modelclass=modelclass)

    @classmethod
    def square(cls, wmin, wmax, transmission):
        """
        Creates a filter with a rectangular shape
        Returns
        -------

        """
        if isinstance(wmin, u.Quantity) is False:
            wmin = wmin*u.AA
        if isinstance(wmax, u.Quantity) is False:
            wmax = wmax*u.AA

        center = (wmax + wmin) * 0.5
        width = wmax - wmin
        modelclass = SpectralElement(Box1D, amplitude=transmission, x_0=center, width=width)

        return cls(modelclass=modelclass)

    def _from_svo(self, filter_name):

        wave, trans = download_svo_filter(filter_name)
        meta = {"filter_name": filter_name, "source": "SVO"}
        self.filter_name = filter_name

        return meta, wave, trans


class ExtinctionCurve(ReddeningLaw, ExtCurveContainer):
    """
    Extinction curves

    TODO: Implement proper __add__, __sub__, __mul__, etc

    """
    def __init__(self, curve_name=None, modelclass=None, **kwargs):

        if curve_name is not None:
            if curve_name in DEFAULT_CURVES:
                curve_name = DEFAULT_CURVES[curve_name]

            ExtCurveContainer.__init__(self, curve_name)
            meta, wave, rvs = self._loader()
            ReddeningLaw.__init__(self, Empirical1D, points=wave, lookup_table=rvs, meta=meta)

        elif modelclass is not None:
            ReddeningLaw.__init__(self, modelclass)
        else:
            raise ValueError("please define a filter")

    @classmethod
    def from_vectors(cls, waves, flux, meta=None, wave_unit=u.AA):
        """
        Create a ``Passband`` directly from from vectos (lists, numpy.arrays, etc)
        Parameters
        ----------
        waves: list-like
        flux: list-like
        meta: dictionary containing the metadata
        wave_unit: u.Quantity, defaulted to angstroms
        flux_unit: u.Quantiy, defaulted to FLAM

        Returns
        -------
        Passband
        """
        if isinstance(waves, u.Quantity) is False:
            waves = waves * wave_unit

        modelclass = SpectralElement(Empirical1D, points=waves, lookup_table=flux, meta=meta)

        return cls(modelclass=modelclass)

    @classmethod
    def from_file(cls, filename, **kwargs):

        modelclass = ReddeningLaw.from_file(filename, **kwargs)

        return cls(modelclass=modelclass)

    def _loader(self):
        """
        Load a filter from the database

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """

        try:  # it should also try to read it from the file directly
            self.wave_unit = units.validate_unit(self.wave_unit)
        except exceptions.SynphotError:
            self.wave_unit = u.AA

        if self.data_type == "fits":
            meta, lam, rvs = read_fits_spec(self.filename, ext=1,
                                            wave_unit=self.wave_unit, flux_unit=self.extinction_unit,
                                            # self._internal_flux_unit,
                                            wave_col=self.wave_column, flux_col=self.extinction_column)
        elif self.data_type == "ascii":
            meta, lam, rvs = read_ascii_spec(self.filename,
                                             wave_unit=self.wave_unit, flux_unit=self.extinction_unit,
                                             # self._internal_flux_unit,
                                             wave_col=self.wave_column, flux_col=self.extinction_column)

        return meta, lam, rvs


class Spextrum(SpectrumContainer, SourceSpectrum):
    """
    Class to handle spectra. This class download, load, stores and manipulates the spectra.

    This class can be initialized with a remote file which will be downloaded from
    the database or with a synphot.BaseSpectrum or SourceSpectrum


    Parameters
    ----------
    template_name : Name of the template to download with format library/template e.g. "kc96/s0"
    modelclass : SourceSpectrum or BaseSpectrum

    """

    def __init__(self, template_name=None, modelclass=None, **kwargs):

        if template_name is not None:
            if template_name in DEFAULT_SPECTRA:
                template_name = DEFAULT_SPECTRA[template_name]

            SpectrumContainer.__init__(self, template_name)
            meta, lam, flux = self._loader()
            SourceSpectrum.__init__(self, Empirical1D, points=lam, lookup_table=flux, meta=meta,   **kwargs)
          #  modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)
        elif modelclass is not None:
           # modelclass = SourceSpectrum(modelclass=modelclass)
            SourceSpectrum.__init__(self, modelclass, **kwargs)
        else:
            raise ValueError("please define a spectra")

#        SourceSpectrum.__init__(self, modelclass, **kwargs)

    def _loader(self):
        """
        Load a template from the database

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """
        try:
            self.wave_unit = units.validate_unit(self.wave_unit)
            self.flux_unit = units.validate_unit(self.flux_unit)
        except exceptions.SynphotError:   # Assumes angtroms and FLAM
            self.wave_unit = u.AA
            self.flux_unit = units.FLAM

     #   self.resolution = self.resolution * self.wave_unit

        if self.data_type == "fits":  # make try and except here to catch most problems
            meta, lam, flux = read_fits_spec(self.path, ext=1,
                                             wave_unit=self.wave_unit, flux_unit=self.flux_unit,
                                             wave_col=self.wave_column_name, flux_col=self.flux_column_name)
        else:
            meta, lam, flux = read_ascii_spec(self.path,
                                              wave_unit=self.wave_unit, flux_unit=self.flux_unit)

        return meta, lam, flux

    @classmethod
    def from_vectors(cls, waves, flux, meta=None, wave_unit=u.AA, flux_unit=units.FLAM):
        """
        Create a ``Passband`` directly from from vectos (lists, numpy.arrays, etc)
        Parameters
        ----------
        waves: list-like
        flux: list-like
        meta: dictionary containing the metadata
        wave_unit: u.Quantity, defaulted to angstroms
        flux_unit: u.Quantiy, defaulted to FLAM

        Returns
        -------
        Passband
        """
        if isinstance(waves, u.Quantity) is False:
            waves = waves * wave_unit
        if isinstance(flux, (u.Quantity, u.core.Unit)) is False:
            flux = flux * flux_unit

        modelclass = SpectralElement(Empirical1D, points=waves, lookup_table=flux, meta=meta)

        return cls(modelclass=modelclass)

    @classmethod
    def from_file(cls, filename, **kwargs):
        """
        This is just an adaptation of the ``synphot.SourceSpectrum.from_file()`` method

        Parameters
        ----------
        filename
        kwargs

        Returns
        -------
        Spextrum instance

        """
        modelclass = SourceSpectrum.from_file(filename, **kwargs)

        return cls(modelclass=modelclass)

    @classmethod
    def from_specutils(cls, spectrum_object):
        """
        This function _tries_ to create a Spectrum from a ``specutils.Spectrum1D`` instance.

        ``specutils.Spectrum1D``  can read multiple file formats with the ``.read`` method.
        please read ``specutils`` documentation.

        Parameters
        ----------
        spectrum_object: specutils.Spectrum1D object

        Returns
        -------
        a Spextrum instance
        """
        meta = spectrum_object.meta
        lam = spectrum_object.spectral_axis
        flux = spectrum_object.flux
        modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)

        return cls(modelclass=modelclass)

    def redshift(self, z=0, vel=0):
        """
        Redshift or blueshift a spectra

        Parameters
        ----------
        z: redshift or
        vel: radial velocity,  if no unit are present it is assumed to be in m/s

        Returns
        -------
        a new Spextrum instance

        """
        if vel != 0:
            if isinstance(vel, u.Quantity) is False:
                vel = vel * u.m / u.s  # assumed to be in m/s

            z = (vel.to(u.m / u.s) / speed_of_light).value
        if z <= -1:
            raise ValueError("Redshift or velocity unphysical")

        lam = self.waveset * (1 + z)
        flux = self(self.waveset)
        self.meta.update({"redshift": z})
        modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=self.meta)
        sp = self._restore_attr(Spextrum(modelclass=modelclass))

        return sp

    def add_emi_lines(self, center, flux, fwhm):
        """
        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)
        Parameters
        ----------
        spectrum: a synphot spectrum
        center: center of the line, astropy.units accepted
        flux: total flux of the line, astropy.units accepted
        fwhm: fwhm of the line, astropy.units accepted

        Returns
        -------
        the spectrum with the emission lines

        """
        if isinstance(center, u.Quantity) is True:
            center = center.to(u.AA).value
        if isinstance(flux, u.Quantity) is True:
            flux = flux.to(u.erg / (u.cm ** 2 * u.s)).value
        if isinstance(fwhm, u.Quantity) is True:
            fwhm = fwhm.to(u.AA).value

        center = np.array([center]).flatten()
        flux = np.array([flux]).flatten()
        fwhm = np.array([fwhm]).flatten()
        sp = self #Spextrum(modelclass=self.model)

        meta = self.meta

        for c, x, f in zip(center, flux, fwhm):
            self.meta.update({"flux_line_" + str(c): str(x)})
            self.meta.update({"fwhm_line_" + str(c): str(f)})

            line = GaussianFlux1D(mean=c, total_flux=x, fwhm=f)
            lam = line.sampleset(factor_step=0.35)  # bit better than Nyquist
            g_em = SourceSpectrum(Empirical1D, points=lam, lookup_table=line(lam))
            sp = sp + g_em

        return sp

    def add_abs_lines(self, center, ew, fwhm):
        """
        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)
        Add a absorption line of to a spectrum with center, fwhm and equivalent width specified by the user
        It also supports emission lines if ew is negative

        Parameters
        ----------
        center
        fwhm
        ew

        Returns
        -------
        a synphot.SourceSpectrum
        """
        if isinstance(center, u.Quantity) is True:
            center = center.to(u.AA).value
        if isinstance(ew, u.Quantity) is True:
            ew = ew.to(u.AA).value
        if isinstance(fwhm, u.Quantity) is True:
            fwhm = fwhm.to(u.AA).value

        center = np.array([center]).flatten()
        ew = np.array([ew]).flatten()
        fwhm = np.array([fwhm]).flatten()

        sp = self  #  .__class__(modelclass=self.model)
        for c, e, f in zip(center, ew, fwhm):
            sign = -1 * np.sign(e)  # to keep the convention that EL are negative and ABS are positive
            left, right = c - np.abs(e / 2), c + np.abs(e / 2)
            wavelengths = self.waveset[(self.waveset.value >= left) & (self.waveset.value <= right)]
            fluxes = units.convert_flux(wavelengths=wavelengths,
                                        fluxes=self(wavelengths),
                                        out_flux_unit=units.FLAM)
            flux = np.trapz(fluxes.value, wavelengths.value)
            line = GaussianFlux1D(mean=c, total_flux=sign * flux, fwhm=f)
            lam = line.sampleset(factor_step=0.35)  # bit better than Nyquist
            g_abs = SourceSpectrum(Empirical1D, points=lam, lookup_table=line(lam))

            sp = sp + g_abs
            self.meta.update({"EW_line_" + str(c): str(e)})
            self.meta.update({"fwhm_line_" + str(c): str(f)})

            if (sp(wavelengths).value < 0).any():
                warnings.warn("Warning: Flux<0 for specified EW and FHWM, setting it to Zero")
                waves = sp.waveset[sp(sp.waveset) < 0]
                zero_sp = SourceSpectrum(Empirical1D, points=waves, lookup_table=-1 * sp(waves).value)
                sp = sp + zero_sp  # Spextrum(modelclass=sp.model + zero_sp.model)

        return sp

    def redden(self, curve_name,  Ebv=0, Av=None, Rv=3.1):
        """
        This function attenuate  a spectrum with a extinction curve normalized to a E(B-V)

        Parameters
        ----------
        curve_name: str, name of the extinction curve
        Av: float, Av parameter
        Ebv: float, E(B-V) color excess

        Returns
        -------
        an attenuated synphot spectrum

        """
        if Av is not None:
            Ebv = Av / Rv

        extcurve = ExtinctionCurve(curve_name)

        extinction = extcurve.extinction_curve(Ebv)
        sp = Spextrum(modelclass=self.model * extinction.model)
        sp = self._restore_attr(sp)
        sp.meta.update({"Reddening with": curve_name, "Reddening with E(B-V)": Ebv})

        return sp

    def deredden(self, curve_name, Av=None, Ebv=0, Rv=3.1):
        """
        This function de-redden a spectrum.

        Parameters
        ----------
        curve_name: str, name of the extinction curve
        Av: float, Av parameter
        Ebv: float, E(B-V) color excess

        Returns
        -------
        an attenuated synphot spectrum

        """
        if Av is not None:
            Ebv = Av / Rv

        ext_curve = ExtinctionCurve(curve_name)

        extinction = ext_curve.extinction_curve(Ebv)
        sp = Spextrum(modelclass=self.model / extinction.model)
        sp = self._restore_attr(sp)
        sp.meta.update({"Dereddening with": curve_name, "Dereddening with E(B-V)": Ebv})

        return sp

    def rebin_spectra(self, new_waves):
        """
        Rebin a synphot spectra to a new wavelength grid conserving flux.
        Grid does not need to be linear and can be at higher or lower resolution

        Parameters
        ----------
        new_waves: an array of the output wavelenghts in Angstroms but other units can be
            specified

        Returns
        -------
        a new Spextrum instance

        """
        if isinstance(new_waves, u.Quantity):
            new_waves = new_waves.to(u.AA).value

        waves = self.waveset.value   # else assumed to be angstroms
        f = np.ones(len(waves))
        filt = SpectralElement(Empirical1D, points=waves, lookup_table=f)
        obs = Observation(self, filt, binset=new_waves, force='taper')
        newflux = obs.binflux
        rebin_spec = Empirical1D(points=new_waves, lookup_table=newflux, meta=self.meta)
        sp = Spextrum(modelclass=rebin_spec)
        sp = self._restore_attr(sp)
        sp.meta.update({"rebinned_spectra": True})

        return sp

    def logrebin(self):
        """
        TODO: Estimate optimal rebinning factors here.
        Returns
        -------
        a Spextrum with a log rebinned in wavelength
        """

        waves = self.waveset.value
        steps = waves[1:] - waves[:-1]
        min_step = np.min(steps)
        wmin, wmax = np.min(waves), np.max(waves)
        print(wmin, wmax, min_step)
        nbins = int((wmax - wmin)/min_step)
        logwaves = np.geomspace(wmin, wmax, nbins)

        return self.rebin_spectra(logwaves)

    def smooth(self, sigma):
        """
        Smooth the Spectrum with a Gaussian Kernel, expressed in km/s (or other velocity unit).

        The spectra is first log-rebinned so it has a constant velocity spacing and after that
        is smoothed with a gaussian kernel


        Parameters
        ----------
        sigma: u.Quantity, width of the Kernel in

        Returns
        -------
        smoothed spextrum
        """
        if isinstance(sigma, u.Quantity) is False:
            sigma = sigma * u.km / u.s

        sigma = sigma.to(u.km / u.s)

        sp_log = self.logrebin()
        lam = sp_log.waveset.value
        flux = sp_log(sp_log.waveset)
        steps = (lam[1:] - lam[:-1]) / lam[:1]
        vel_steps = steps * speed_of_light.to(u.km/u.s)
        step_size = np.median(vel_steps)

        if step_size > sigma:
            warnings.warn("spectra is undersampled for the provided sigma value")

        conv_sigma = sigma / step_size
        smoothed_flux = gaussian_filter1d(flux, conv_sigma.value)

        self.meta.update({"KERNEL_SIZE": sigma.value})
        modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=smoothed_flux, meta=self.meta)
        sp = Spextrum(modelclass=modelclass)
        sp = sp._restore_attr(sp)

        return sp

    @classmethod
    def flat_spectrum(cls, mag=0, system_name="AB", wavelengths=None):
        """
        Creates a flat spectrum in the preferred system scaled to a magnitude,
        default a zero magnitude spectrum
        Parameters
        ----------
        mag: float,
            magnitude of the reference spectrum, default=0
        system_name: AB, Vega or ST, default AB

        wavelengths: The waveset of the reference spectrum if not Vega

        Returns
        -------
        a Spextrum instance
        """
        if wavelengths is None:  # set a default waveset with R~805
            wavelengths, info = utils.generate_wavelengths(minwave=100, maxwave=50000, num=5000,
                                                           log=True, wave_unit=u.AA)
        if system_name.lower() in ["vega"]:
            spec = get_vega_spectrum()
            spec = spec * 10**(-0.4*mag)
        elif system_name.lower() in ["ab"]:
            spec = SourceSpectrum(ConstFlux1D, amplitude=mag * u.ABmag)
            spec = SourceSpectrum(Empirical1D, points=wavelengths,
                                  lookup_table=spec(wavelengths, flux_unit=u.ABmag))
        elif system_name.lower() in ["st", "hst"]:
            spec = SourceSpectrum(ConstFlux1D, amplitude=mag * u.STmag)
            spec = SourceSpectrum(Empirical1D, points=wavelengths,
                                  lookup_table=spec(wavelengths, flux_unit=u.STmag))
        else:
            raise ValueError("only AB, ST and Vega systems are supported")

        return cls(modelclass=spec)

    @classmethod
    def black_body_spectrum(cls, temperature=9500, amplitude=0, filter_name=None, filter_file=None):
        """
        Produce a blackbody spectrum for a given temperature and scale it to a magnitude
        in a filter

        Parameters
        ----------
        temperature: the temperature in Kelvin degrees
        amplitude: `astropy.Quantity``, float
                The value that the spectrum should have in the given filter. Acceptable
                astropy quantities are:
                - u.mag : Vega magnitudes
                - u.ABmag : AB magnitudes
                - u.STmag : HST magnitudes
                - u.Jy : Jansky per filter bandpass
                Additionally the ``FLAM`` and ``FNU`` units from ``synphot.units`` can
                be used when passing the quantity for ``amplitude``:

        filter_name : str
                Name of a filter from
                - a generic filter name (see ``FILTER_DEFAULTS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database

        filter_file: str
                A file with a transmission curve

        Returns
        -------
        a scaled black-body spectrum
        """
        sp = cls(modelclass=BlackBody1D, temperature=temperature)

        return sp.scale_to_magnitude(amplitude=amplitude,
                                     filter_name=filter_name,
                                     filter_file=filter_file)

    @classmethod
    def powerlaw(cls, alpha=1, amplitude=0, filter_name=None, filter_file=None):
        """
        Return a power law spectrum F(lambda) ~ lambda^alpha scaled to a magnitude
        (amplitude) in an particular band

        Parameters
        ----------
        alpha
        amplitude
        filter_name
        filter_file

        TODO: definition of x0 as filter pivot wavelength
        Returns
        -------

        """

        sp = cls(modelclass=PowerLawFlux1D, amplitude=amplitude, x_0=2000, alpha=alpha)

        return sp.scale_to_magnitude(amplitude=amplitude, filter_name=filter_name, filter_file=filter_file)

    def scale_to_magnitude(self, amplitude, filter_name=None, filter_file=None):
        """
            Scales a Spectrum to a value in a filter
            copied from scopesim.effects.ter_curves.scale_spectrum with slight modifications
            Parameters
            ----------
            amplitude : ``astropy.Quantity``, float
                The value that the spectrum should have in the given filter. Acceptable
                astropy quantities are:
                - u.mag : Vega magnitudes
                - u.ABmag : AB magnitudes
                - u.STmag : HST magnitudes
                - u.Jy : Jansky per filter bandpass
                Additionally the ``FLAM`` and ``FNU`` units from ``synphot.units`` can
                be used when passing the quantity for ``amplitude``:

            filter_name : str
                Name of a filter from
                - a generic filter name (see ``FILTER_DEFAULTS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
            filter_file: str
                A file with a transmission curve


            Returns
            -------
            spectrum : a Spectrum
                Input spectrum scaled to the given amplitude in the given filter
        """

        if filter_file is not None:
            filter_curve = Passband.from_file(filename=filter_file)
        else:
            filter_curve = Passband(filter_name=filter_name)

        if isinstance(amplitude, u.Quantity):
            if amplitude.unit.physical_type == "spectral flux density":
                if amplitude.unit != u.ABmag:
                    amplitude = amplitude.to(u.ABmag)
                ref_spec = self.flat_spectrum(mag=amplitude.value, system_name="AB")

            elif amplitude.unit.physical_type == "spectral flux density wav":
                if amplitude.unit != u.STmag:
                    amplitude = amplitude.to(u.STmag)
                ref_spec = self.flat_spectrum(mag=amplitude.value, system_name="ST")

            elif amplitude.unit == u.mag:
                ref_spec = self.flat_spectrum(mag=amplitude.value, system_name="Vega")

            else:
                raise ValueError("Units of amplitude must be one of "
                                 "[u.mag, u.ABmag, u.STmag]: {}".format(amplitude))
        else:
            ref_spec = self.flat_spectrum(mag=amplitude, system_name="Vega")

        ref_flux = Observation(SourceSpectrum(modelclass=ref_spec),
                               filter_curve).effstim(flux_unit=units.PHOTLAM)
        real_flux = Observation(SourceSpectrum(modelclass=self),
                                filter_curve).effstim(flux_unit=units.PHOTLAM)

        scale_factor = ref_flux / real_flux
        sp = self * scale_factor

        return sp

    def get_magnitude(self, filter_name=None, filter_file=None, system_name="AB"):
        """
            Obtain the magnitude in filter for a user specified photometric system

            Parameters
            ----------
            filter_name : str
                Name of a filter from
                - a generic filter name (see ``FILTER_DEFAULTS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
            filter_file: str
                A file with a transmission curve
            system_name: str
                The photometric system Vega, AB or ST


            Returns
            -------
            spectrum : a Spextrum
                Input spectrum scaled to the given amplitude in the given filter
        """

        if filter_file is not None:
            filter_curve = Passband.from_file(filename=filter_file)
        else:
            filter_curve = Passband(filter_name=filter_name)

        if system_name.lower() in ["vega"]:
            unit = u.mag
        elif system_name.lower() in ["st", "hst"]:
            unit = u.STmag
        else:
            unit = u.ABmag

        ref_spec = self.flat_spectrum(mag=0, system_name=system_name)
        ref_flux = Observation(SourceSpectrum(modelclass=ref_spec),
                               filter_curve).effstim(flux_unit=units.PHOTLAM)
        real_flux = Observation(self, filter_curve).effstim(flux_unit=units.PHOTLAM)
        mag = -2.5*np.log10(real_flux.value/ref_flux.value)

        return mag * unit

    def photons_in_range(self, wmin=None, wmax=None, area=1*u.cm**2,
                         filter_name=None, filter_file=None):
        """
        Return the number of photons between wave_min and wave_max or within
        a bandpass (filter)

        Parameters
        ----------
        wmin :
            [Angstrom]
        wmax :
            [Angstrom]
        area : u.Quantity
            [cm2]
        filter_name :
        filter_file :

        Returns
        -------
        counts : u.Quantity array

        """
        if isinstance(area, u.Quantity):
            area = area.to(u.cm ** 2).value  #
        if isinstance(wmin, u.Quantity):
            wmin = wmin.to(u.Angstrom).value
        if isinstance(wmax, u.Quantity):
            wmin = wmax.to(u.Angstrom).value

        if (filter_name is None) and (filter_file is None):
            # this makes a bandpass out of wmin and wmax
            try:
                mid_point = 0.5 * (wmin + wmax)
                width = abs(wmax - wmin)
                filter_curve = SpectralElement(Box1D, amplitude=1, x_0=mid_point, width=width)
            except ValueError("Please specify wmin/wmax or a filter"):
                raise

        elif filter_file is not None:
            filter_curve = Passband.from_file(filename=filter_file)
        else:
            filter_curve = Passband(filter_name=filter_name)

        obs = Observation(self, filter_curve)
        counts = obs.countrate(area=area * u.cm ** 2)

        return counts

    def add_noise(self, wmin, wmax):
        """
        Returns a spectra with a given S/N in the specified wavelength range
        Returns
        -------

        """
        pass

    def _restore_attr(self, spextrum):
        """
        just a hack to down stream important attributes in particular methods.
        there must be a better way
        """
        temp_dict = {k: self.__dict__[k] for k in self.__dict__ if k.startswith("_") is False}
        spextrum.__dict__.update(temp_dict)
        return spextrum


    # ------ Copied from synphot.SourceSpectrum so operations can also happen here --------

    def __add__(self, other):
        """Add ``self`` with ``other``."""
        self._validate_other_add_sub(other)
        result = self.__class__(modelclass=self.model + other.model)
        self._merge_meta(self, other, result)
        return result

    def __sub__(self, other):
        """Subtract other from self."""
        self._validate_other_add_sub(other)
        result = self.__class__(modelclass=self.model - other.model)
        self._merge_meta(self, other, result)
        return result

    def __mul__(self, other):
        """Multiply self and other."""
        self._validate_other_mul_div(other)

        if isinstance(other, (u.Quantity, numbers.Number)):
            newcls = self.__class__(modelclass=self.model | Scale(other))
        elif isinstance(other, BaseUnitlessSpectrum):
            newcls = self.__class__(modelclass=self.model * other.model)
        else:  # Source spectrum
            raise exceptions.IncompatibleSources(
                'Cannot multiply two source spectra together')

        self._merge_meta(self, other, newcls)
        return newcls

    def _validate_other_add_sub(self, other):
        """
        Conditions for other to satisfy before add/sub.
        Copied from synphot.SourceSpectrum so additions can happen with SourceSpectrum
        """
        if not isinstance(other, (self.__class__, SourceSpectrum)):
            raise exceptions.IncompatibleSources(
                'Can only operate on {0}.'.format(self.__class__.__name__))

    def _validate_other_mul_div(self, other):
        """Conditions for other to satisfy before mul/div.
           Copied from synphot.SourceSpectrum so additions can happen with a SourceSpectrum

        """

        if not isinstance(other, (u.Quantity, numbers.Number,
                                  BaseUnitlessSpectrum, SourceSpectrum, self.__class__)):
            raise exceptions.IncompatibleSources(
                'Can only operate on scalar number/Quantity or spectrum')
        elif (isinstance(other, u.Quantity) and
              (other.unit.decompose() != u.dimensionless_unscaled or
               not np.isscalar(other.value) or
               not isinstance(other.value, numbers.Real))):
            raise exceptions.IncompatibleSources(
                'Can only operate on real scalar dimensionless Quantity')
        elif (isinstance(other, numbers.Number) and
              not (np.isscalar(other) and isinstance(other, numbers.Real))):
            raise exceptions.IncompatibleSources(
                'Can only operate on real scalar number')


#------------------------------ END    -------------------------------------------


def get_vega_spectrum():
    """
    Retrieve the Vega spectrum from the database

    Notes
    -----
    To access wavelength and fluxes use::
        sp = get_vega_spectrum()
        waves, fluxes = sp.waveset, sp(sp.waveset)

    """
    return Spextrum("ref/vega")




