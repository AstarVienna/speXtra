# -*- coding: utf-8 -*-
"""
speXtra: A python tool to manage and manipulate astronomical spectra
"""
import numbers
import warnings
import os

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
    Class to handle astronomical filters

    TODO: Implement proper __add__, __sub__, __mul__, etc

    """

    def __init__(self, filter_name=None, modelclass=None, **kwargs):

        if filter_name is not None:
            if filter_name in DEFAULT_FILTERS:
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
            SpectralElement.__init__(self, modelclass, **kwargs)
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
            meta, lam, trans = read_fits_spec(self.filename, ext=1, wave_unit=self.wave_unit)
                                        #      wave_col=self.wave_column_name, flux_col=self.trans_column_name)
        elif self.data_type == "ascii":
            meta, lam, trans = read_ascii_spec(self.filename,  wave_unit=self.wave_unit)
                                     #          wave_col=self.wave_column_name, flux_col=self.trans_column_name)
                                             # , flux_unit=self._internal_flux_unit,
        lam = lam.to(u.AA)

        return meta, lam, trans.value

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
        self.wave_unit = u.AA

        return meta, wave, trans


class ExtinctionCurve(ReddeningLaw, ExtCurveContainer):
    """
    Class to handle extinction curves

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
            SourceSpectrum.__init__(self, Empirical1D, points=lam, lookup_table=flux, meta=meta,  **kwargs)

            self.repr = "Spextrum(%s)" % self.template
        elif modelclass is not None:

            SourceSpectrum.__init__(self, modelclass=modelclass, **kwargs)

        else:
            raise ValueError("please define a spectra")


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
    def from_arrays(cls, waves, flux, meta=None, wave_unit=u.AA, flux_unit=units.FLAM):
        """
        Create a ``Passband`` directly from from arrays (lists, numpy.arrays, etc)
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

        modelclass = SourceSpectrum(Empirical1D, points=waves, lookup_table=flux, meta=meta)

        sp = cls(modelclass=modelclass)
        sp.repr = repr(sp.model)
        return sp

    @classmethod
    def from_file(cls, filename, **kwargs):
        """
        This is just an wrapper of the ``synphot.SourceSpectrum.from_file()`` method

        Parameters
        ----------
        filename
        kwargs

        Returns
        -------
        Spextrum instance

        """
        modelclass = SourceSpectrum.from_file(filename, **kwargs)
        sp = cls(modelclass=modelclass)
        sp.repr = "Spextrum.from_file(%s)" % filename

        return sp

    @classmethod
    def flat_spectrum(cls, amplitude=0, waves=None):
        """
        Creates a flat spectrum in the preferred system scaled to a magnitude,
        default a zero magnitude spectrum
        Parameters
        ----------
        amplitude: float, u.Quantity
            amplitude/magnitude of the reference spectrum, default=0
            default is u.ABmag
            for vega use u.mag or ``synphot.units.VEGAMAG``

        waves: The waveset of the reference spectrum if not Vega
           if not provided they will be created at a resolution of R~800

        Returns
        -------
        a Spextrum instance
        """
        if waves is None:  # set a default waveset with R~805
            waves, info = utils.generate_wavelengths(minwave=100, maxwave=50000, num=5000,
                                                     log=True, wave_unit=u.AA)

        if isinstance(amplitude, u.Quantity) is False:
            amplitude = amplitude*u.ABmag

        if amplitude.unit is u.Unit("mag") or amplitude.unit is u.Unit("vegamag"):
            spec = get_vega_spectrum()
            spec = spec * 10 ** (-0.4 * amplitude.value)
            system_name = amplitude.unit
        else:
            const = ConstFlux1D(amplitude=amplitude)
            spec = cls(modelclass=Empirical1D, points=waves.value, lookup_table=const(waves.value))
            system_name = amplitude.unit

        spec.repr = "Spextrum.flat_spectrum(amplitude=%s)" % str(amplitude)

        return spec

    @classmethod
    def black_body_spectrum(cls, temperature=9500, amplitude=0, filter_curve=None, waves=None):
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

        filter_curve : str
                Name of a filter from
                - a generic filter name (see ``FILTER_DEFAULTS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
                - a filename with the filter file
                - a ``Passband`` or ``synphot.SpectralElement`` object


        Returns
        -------
        a scaled black-body spectrum
        """

        if waves is None:  # set a default waveset with R~805
            waves, info = utils.generate_wavelengths(minwave=100, maxwave=50000, num=5000,
                                                     log=True, wave_unit=u.AA)

        if isinstance(amplitude, u.Quantity) is False:
            amplitude = amplitude*u.ABmag

        bb = BlackBody1D(temperature=temperature)

        sp = cls(modelclass=Empirical1D, points=waves, lookup_table=bb(waves))
        sp = sp.scale_to_magnitude(amplitude=amplitude, filter_curve=filter_curve)
        sp.repr = "Spextrum.black_body_spectrum(amplitude=%s, temperature=%s, filter_curve=%s)" % \
                  (str(amplitude), str(temperature), str(filter_curve))

        return sp

    @classmethod
    def powerlaw(cls, alpha=1, x_0=5000, amplitude=0, filter_curve=None, waves=None):
        """
        Return a power law spectrum F(lambda) ~ lambda^alpha scaled to a magnitude
        (amplitude) in an particular band

        Parameters
        ----------
        alpha : The spectral slope
        x_0 : float, u.Quantity
           Pivot wavelength
        amplitude : float, u.Quantity
            normalize the spectrum to that quantity
        filter_curve : str
                Name of a filter from
                - a generic filter name (see ``FILTER_DEFAULTS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
                - a filename with the filter file
                - a ``Passband`` or ``synphot.SpectralElement`` object

        Returns
        -------

        """
        if waves is None:  # set a default waveset with R~805
            waves, info = utils.generate_wavelengths(minwave=100, maxwave=50000, num=5000,
                                                     log=True, wave_unit=u.AA)

        if isinstance(amplitude, u.Quantity) is False:
            amplitude = amplitude*u.ABmag

        pl = SourceSpectrum(PowerLawFlux1D, amplitude=1, x_0=x_0, alpha=alpha)
        sp = cls(modelclass=Empirical1D, points=waves, lookup_table=pl(waves))

        sp = sp.scale_to_magnitude(amplitude=amplitude, filter_curve=filter_curve)

        sp.repr = "Spextrum.powerlaw(alpha=%s, x_0=%s, amplitude=%s, filter_curve=%s)" % \
                  (str(alpha), str(x_0), str(amplitude), str(filter_curve))

        return sp

    @classmethod
    def emission_line_spectra(cls, center, fwhm,  flux,  amplitude=40*u.ABmag, waves=None):
        """
        Create a emission line spextrum superimpossed to a faint continuum


        Parameters
        ----------
        center
        fwhm
        flux
        amplitude
        waves

        Returns
        -------

        """
        if isinstance(center, u.Quantity) is False:
            center = center * u.AA
        else:
            center = center.to(u.AA, equivalencies=u.spectral())
        if isinstance(fwhm, u.Quantity) is False:
            fwhm = fwhm * u.AA
        else:
            fwhm = fwhm.to(u.AA, equivalencies=u.spectral())

        if waves is None:
            wmin = np.min(center.value) - 2000
            wmax = np.max(center.value) + 2000
            step = np.min(fwhm.value) / 3   # for Nyquist sampling
            waves = np.arange(wmin, wmax, step) * u.AA

        sp = cls.flat_spectrum(amplitude=amplitude, waves=waves)
        sp.add_emi_lines(center=center, fwhm=fwhm, flux=flux)

        return sp

    @property
    def wave_min(self):
        return np.min(self.waveset)

    @property
    def wave_max(self):
        return np.max(self.waveset)

    def cut(self, wave_min, wave_max):
        """
        Cut the spectrum between wmin and wmax

        Parameters
        ----------
        wave_min: float, u.Quantity,
        wave_max: float, u.Quantity,

        Returns
        -------

        a new spextrum
        """
        if isinstance(wave_min, u.Quantity) is False:
            wave_min = wave_min * u.AA
        if isinstance(wave_max, u.Quantity) is False:
            wave_max = wave_max * u.AA

        wave_min = wave_min.to(u.AA)
        wave_max = wave_max.to(u.AA)
        new_waves = self.waveset[(self.waveset >= wave_min) & (self.waveset <= wave_max)]
        sp = Spextrum(modelclass=Empirical1D, points=new_waves, lookup_table=self(new_waves), meta=self.meta)
        sp = self._restore_attr(sp)

        return sp

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
        sp = Spextrum(modelclass=Empirical1D, points=lam, lookup_table=flux, meta=self.meta)
        sp = self._restore_attr(sp)

        return sp

    def scale_to_magnitude(self, amplitude, filter_curve=None):
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

            filter_curve : str
                Name of a filter from
                - a generic filter name (see ``DEFAULT_FILTERS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
                - the path to a filename containing the curve (see ``Passband``)
                - a ``Passband`` or ``synphot.SpectralElement`` instance
            filter_file: str
                A file with a transmission curve


            Returns
            -------
            spectrum : a Spectrum
                Input spectrum scaled to the given amplitude in the given filter
        """

        if os.path.exists(filter_curve):
            filter_curve = Passband.from_file(filename=filter_curve)
        elif isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = filter_curve
        else:
            filter_curve = Passband(filter_curve)

        if isinstance(amplitude, u.Quantity) is False:
            amplitude = amplitude*u.ABmag


       #     if amplitude.unit.physical_type == "spectral flux density":
       #         if amplitude.unit != u.ABmag:
       #             amplitude = amplitude.to(u.ABmag)
       #         ref_spec = self.flat_spectrum(amplitude=amplitude.value, system_name="AB")

        #    elif amplitude.unit.physical_type == "spectral flux density wav":
       #         if amplitude.unit != u.STmag:
       #             amplitude = amplitude.to(u.STmag)
       #         ref_spec = self.flat_spectrum(amplitude=amplitude.value, system_name="ST")

#            elif amplitude.unit == u.mag:
 #               ref_spec = self.flat_spectrum(amplitude=amplitude.value, system_name="Vega")

  #          else:
   #             raise ValueError("Units of amplitude must be one of "
    #                             "[u.mag, u.ABmag, u.STmag]: {}".format(amplitude))
    #    else:

        ref_spec = self.flat_spectrum(amplitude=amplitude)

        ref_flux = Observation(ref_spec, filter_curve).effstim(flux_unit=units.PHOTLAM)
        real_flux = Observation(self, filter_curve).effstim(flux_unit=units.PHOTLAM)

        scale_factor = ref_flux / real_flux
        sp = self * scale_factor
        sp = self._restore_attr(sp)
        sp.meta.update({"magnitude": amplitude,
                        "filter_curve": filter_curve})

        return sp

    def get_magnitude(self, filter_curve=None,  system_name="AB"):
        """
            Obtain the magnitude in filter for a user specified photometric system

            Parameters
            ----------
            filter_curve : str
                Name of a filter from
                - a generic filter name (see ``DEFAULT_FILTERS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
                - the path to the file containing the filter (see ``Passband``)
                - a ``Passband`` or ``synphot.SpectralElement`` object
            system_name: str
                The photometric system Vega, AB or ST


            Returns
            -------
            spectrum : a Spextrum
                Input spectrum scaled to the given amplitude in the given filter
        """
        if os.path.exists(filter_curve):
            filter_curve = Passband.from_file(filename=filter_curve)
        elif isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = filter_curve
        else:
            filter_curve = Passband(filter_curve)

        if system_name.lower() in ["vega"]:
            unit = u.mag
        elif system_name.lower() in ["st", "hst"]:
            unit = u.STmag
        else:
            unit = u.ABmag

        ref_spec = self.flat_spectrum(amplitude=0*unit)
        ref_flux = Observation(ref_spec, filter_curve).effstim(flux_unit=units.PHOTLAM)
        real_flux = Observation(self, filter_curve).effstim(flux_unit=units.PHOTLAM)

        mag = -2.5*np.log10(real_flux.value/ref_flux.value)

        return mag * unit

    def photons_in_range(self, wmin=None, wmax=None, area=1*u.cm**2,
                         filter_curve=None):
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
        filter_curve : str
                Name of a filter from
                - a generic filter name (see ``DEFAULT_FILTERS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
                - the path to the file containing the filter (see ``Passband``)
                - a ``Passband`` or ``synphot.SpectralElement`` object

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

        if filter_curve is None:
            # this makes a bandpass out of wmin and wmax
            try:
                mid_point = 0.5 * (wmin + wmax)
                width = abs(wmax - wmin)
                filter_curve = SpectralElement(Box1D, amplitude=1, x_0=mid_point, width=width)
            except ValueError("Please specify wmin/wmax or a filter"):
                raise

        elif os.path.exists(filter_curve):
            filter_curve = Passband.from_file(filename=filter_curve)
        elif isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = filter_curve
        else:
            filter_curve = Passband(filter_curve)

        obs = Observation(self, filter_curve)
        counts = obs.countrate(area=area * u.cm ** 2)

        return counts

    def get_flux(self, wmin=None, wmax=None, filter_curve=None, flux_unit=units.FLAM):
        """
        Return the flux within a passband

        Parameters
        ----------
        wmin : float, u.Quantity
           minimum wavelength
        wmax : float, u.Quantity
          maximum wavelength
        filter_curve : str
                Name of a filter from
                - a generic filter name (see ``DEFAULT_FILTERS``)
                - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
                - a filter in the spextra database
                - the path to the file containing the filter (see ``Passband``)
                - a ``Passband`` or ``synphot.SpectralElement`` object
        flux_unit: synphot.units, u.Quantity

        Returns
        -------

        """

        if isinstance(wmin, u.Quantity):
            wmin = wmin.to(u.Angstrom).value
        if isinstance(wmax, u.Quantity):
            wmin = wmax.to(u.Angstrom).value

        if filter_curve is None:
            # this makes a bandpass out of wmin and wmax
            try:
                filter_curve = Passband.square(wmin=wmin, wmax=wmax)

            except ValueError("Please specify wmin/wmax or a filter"):
                raise
        elif os.path.exists(filter_curve):
            filter_curve = Passband.from_file(filename=filter_curve)
        elif isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = filter_curve
        else:
            filter_curve = Passband(filter_curve)

        flux = Observation(self, filter_curve).effstim(flux_unit=flux_unit)

        return flux

    def add_emi_lines(self, center, fwhm, flux):
        """
        Add emission lines to an `Spextrum`

        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)\

        Parameters
        ----------
        center: float, list, np.ndarray, u.Quantity
            The center of the line

        fwhm:  float, list, np.ndarray, u.Quantity
            The FWHM of the line
        flux: float, list, np.ndarray, u.Quantity
            The Equivalent Width of the line

        Returns
        -------
        Spextrum

        """
        if isinstance(center, u.Quantity) is True:
            center = center.to(u.AA, equivalencies=u.spectral()).value
        if isinstance(flux, u.Quantity) is False:
            flux = flux*u.erg / (u.cm ** 2 * u.s)
#            flux = flux.to(u.erg / (u.cm ** 2 * u.s), equivalencies=u.spectral()).value
        if isinstance(fwhm, u.Quantity) is True:
            fwhm = fwhm.to(u.AA, equivalencies=u.spectral()).value

        centers = np.array([center]).flatten()
        fluxes = np.array([flux.value]).flatten()*flux.unit
        fwhms = np.array([fwhm]).flatten()
        sp = self
        sp.meta.update({"em_lines": {"center": list(centers),
                                     "flux": list(fluxes),
                                     "fwhm": list(fwhms)}})

        for c, f, w in zip(centers, fluxes, fwhms):

            line = GaussianFlux1D(mean=c, total_flux=f, fwhm=w)
            lam = line.sampleset(factor_step=0.3)  # bit better than Nyquist
            g_em = SourceSpectrum(Empirical1D, points=lam, lookup_table=line(lam)) #, meta=sp.meta)

            sp = sp + g_em

        sp = self._restore_attr(Spextrum(modelclass=sp))

        return sp

    def add_abs_lines(self, center, ew, fwhm):
        """
        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)
        Add a absorption line of to a spectrum with center, fwhm and equivalent width specified by the user
        It also supports emission lines if ew is negative

        Parameters
        ----------
        center: float, list, np.ndarray, u.Quantity
            The center of the line

        fwhm:  float, list, np.ndarray, u.Quantity
            The FWHM of the line
        ew: float, list, np.ndarray, u.Quantity
            The Equivalent Width of the line

        Returns
        -------
        Spextrum
        """
        if isinstance(center, u.Quantity) is True:
            center = center.to(u.AA).value
        if isinstance(ew, u.Quantity) is True:
            ew = ew.to(u.AA).value
        if isinstance(fwhm, u.Quantity) is True:
            fwhm = fwhm.to(u.AA).value

        centers = np.array([center]).flatten()
        ews = np.array([ew]).flatten()
        fwhms = np.array([fwhm]).flatten()
        sp = self  #  .__class__(modelclass=self.model)

        sp.meta.update({"em_lines": {"center": list(centers),
                                     "ew": list(ews),
                                     "fwhm": list(fwhms)}})

        for c, e, f in zip(centers, ews, fwhms):
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

            if (sp(wavelengths).value < 0).any():
                warnings.warn("Warning: Flux<0 for specified EW and FHWM, setting it to Zero")
                waves = sp.waveset[sp(sp.waveset) < 0]
                zero_sp = SourceSpectrum(Empirical1D, points=waves, lookup_table=-1 * sp(waves).value)
                sp = sp + zero_sp  # Spextrum(modelclass=sp.model + zero_sp.model)

        sp = self._restore_attr(Spextrum(modelclass=sp))

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
        smoothed = Empirical1D(points=lam, lookup_table=smoothed_flux, meta=self.meta)
        sp = self._restore_attr(Spextrum(modelclass=smoothed))

        return sp


    def add_noise(self, wmin, wmax):
        """
        Returns a spectra with a given S/N in the specified wavelength range
        Returns
        -------

        """
        raise NotImplementedError

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

    def __repr__(self):

        rep = "<%s>" % self.repr

        return rep
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




