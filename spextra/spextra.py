# -*- coding: utf-8 -*-
"""
speXtra: A python tool to manage and manipulate astronomical spectra
"""
import numbers
import warnings

import numpy as np
from scipy.ndimage import gaussian_filter1d

from astropy.table import Table
import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.modeling.models import Scale

import synphot
from synphot import (units, SourceSpectrum, SpectralElement, Observation, BaseUnitlessSpectrum)
from synphot.models import (Empirical1D, GaussianFlux1D, Box1D, ConstFlux1D, BlackBody1D)
from synphot.specio import read_ascii_spec, read_fits_spec, read_spec
from synphot import exceptions

from .database import SpectralTemplate, Filter, ExtinctionCurve


__all__ = ["Spextrum", "make_passband",  "get_vega_spectrum"]


def make_passband(filter_name=None, filter_file=None, wave_unit=u.Angstrom):
    """
    Make a SpectralElement (aka a synphot passband) from an
    user specified filter in the database.
    Optionally, specify a file in disk and the wavelength units.

    Parameters
    ----------
    filter_name : str
                  a filter name expressed as ``instrument/filter_name``
    filter_file : str,
                  Optionally, make a pasband from a local file
    wave_unit : u.Quantity, optional
                if not specified in the file
                default: u.Angstrom

    Returns
    -------
    passband : a ``synphot.SpectralElement``
    """
    if filter_file is not None:
        try:
            meta, wave, trans = read_spec(filter_file, wave_unit=wave_unit,
                                          flux_unit='transmission')
        except (FileNotFoundError, exceptions.SynphotError) as e:
            print("File not found or malformed", e)
            raise
    else:
        filt = Filter(filter_name=filter_name)
        path, meta = filt.path, filt.meta
        if meta is None:  # it's a svo filter!
            trans_table = Table.read(path, format="votable")
            wave = trans_table['Wavelength'].data.data * u.Angstrom
            trans = trans_table['Transmission'].data.data
        else:   # filter in spextra database
            wave_unit = units.validate_unit(meta["wave_unit"])
            data_type = meta["data_type"]
            if data_type == "fits":
                trans_table = Table.read(path, format="fits")
                wave = trans_table[0][:].data * wave_unit
                trans = trans_table[1][:].data
            elif data_type == "ascii":
                trans_table = Table.read(path, format="ascii")
                wave = trans_table[0][:].data * wave_unit
                trans = trans_table[1][:].data

    passband = SpectralElement(Empirical1D, points=wave, lookup_table=trans, meta=meta)

    return passband


class Spextrum(SourceSpectrum):
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

        self.template_name = template_name
        self.resolution = None
        self.wmin = None
        self.wmax = None
        self.path = None
        self.data_type = None

        if self.template_name is not None:
            meta, lam, flux = self._loader
            modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)
        if modelclass is not None:
            modelclass = modelclass
        else:
            raise ValueError("please define a spectra")

        super().__init__(modelclass, **kwargs)

    @property
    def _loader(self):
        """
        Load a template from the database

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """
        template = SpectralTemplate(self.template_name)
        location, meta = template.path, template.meta
        self.path = location
        self.data_type = meta["data_type"]

        try:  # it should also try to read it from the file directly
            wave_unit = units.validate_unit(meta["wave_unit"])
            flux_unit = units.validate_unit(meta["flux_unit"])
        except exceptions.SynphotError:
            wave_unit = u.AA
            flux_unit = units.FLAM

        self.resolution = meta["resolution"] * wave_unit
        wave_column_name = meta["wave_column_name"]  # same here
        flux_column_name = meta["flux_column_name"]
        file_extension = meta["file_extension"]

        # make try and except here to catch most problems
        if self.data_type == "fits":
            meta, lam, flux = read_fits_spec(location, ext=1,
                                             wave_unit=wave_unit, flux_unit=flux_unit,
                                             wave_col=wave_column_name, flux_col=flux_column_name)
        else:
            meta, lam, flux = read_ascii_spec(location, wave_unit=wave_unit,
                                              flux_unit=flux_unit)

        return meta, lam, flux

    def spectral_edges(self):
        self.wmin = np.min(self.waveset)
        self.wmax = np.max(self.waveset)

        return self.wmin, self.wmax

    @classmethod
    def from1dspec(cls, filename, format="wcs1d-fits", **kwargs):
        """
        This function _tries_ to create a Spectrum from 1d fits files.
        It relies in specutils for its job
        Parameters
        ----------
        filename: str a filename with the spectra
        format: format of the spectra accepted by `specutils` (see `specutils` documentation)

        Returns
        -------
        a Spextrum instance

        """
        try:
            from specutils import Spectrum1D
        except ImportError as ie:
            print(ie, "specutils not installed, cannot import spectra")
            raise

        spec1d = Spectrum1D.read(filename, format=format, **kwargs)
        meta = spec1d.meta
        lam = spec1d.spectral_axis
        flux = spec1d.flux
        modelclass = SourceSpectrum(Empirical1D,
                                    points=lam, lookup_table=flux, meta=meta)

        return cls(modelclass=modelclass)

    def redshift(self, z=0, vel=0):
        """
        Redshift or blueshift a spectra

        Parameters
        ----------
        z: redshift
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
        meta = self.meta.update({"redshift": z})
        modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)

        return Spextrum(modelclass=modelclass)

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

        sp = self  #.__class__(modelclass=self.model)
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
        curve: the extinction curve
        EBV: E(B-V)

        Returns
        -------
        an attenuated synphot spectrum

        """
        if Av is not None:
            Ebv = Av / Rv

        extcurve = ExtinctionCurve(curve_name)
        curve, meta = extcurve.path, extcurve.meta

        if meta["data_type"] == "fits":
            header, wavelengths, rvs = read_fits_spec(curve, flux_col='Av/E(B-V)')
        else:
            header, wavelengths, rvs = read_ascii_spec(curve)

        red_law = synphot.ReddeningLaw(Empirical1D, points=wavelengths, lookup_table=rvs, meta={'header': header})
        extinction = red_law.extinction_curve(Ebv)
        sp = Spextrum(modelclass=self.model * extinction.model)

        return sp

    def deredden(self, curve_name, Av=None, Ebv=0, Rv=3.1):
        """
        This function de-redden a spectrum.

        Parameters
        ----------
        spectrum: a synphot spectrum
        curve: the extinction curve
        EBV: E(B-V)

        Returns
        -------
        an attenuated synphot spectrum

        """
        if Av is not None:
            Ebv = Av / Rv

        ext_curve = ExtinctionCurve(curve_name)
        curve, meta = ext_curve.path, ext_curve.meta

        if meta["data_type"] == "fits":
            header, wavelengths, rvs = read_fits_spec(curve, flux_col='Av/E(B-V)')
        else:
            header, wavelengths, rvs = read_ascii_spec(curve)

        red_law = synphot.ReddeningLaw(Empirical1D, points=wavelengths, lookup_table=rvs, meta={'header': header})
        extinction = red_law.extinction_curve(Ebv)
        sp = Spextrum(modelclass=self.model / extinction.model)

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

        return Spextrum(modelclass=rebin_spec)

    def logrebin(self):
        """
        TODO: Estimate optimal rebinning factors here.
        Returns
        -------
        a Spectrum with a log rebinned in wavelength
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

        flux_unit = flux.unit
        meta = self.meta

        meta.update({"KERNEL_SIZE": sigma.value})
        modelclass = SourceSpectrum(Empirical1D,
                                    points=lam, lookup_table=smoothed_flux,
                                    meta=meta)

        return Spextrum(modelclass=modelclass)

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
            wavelengths, info = synphot.utils.generate_wavelengths(minwave=100, maxwave=50000, num=5000,
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
            filter_curve = make_passband(filter_file=filter_file)
        else:
            filter_curve = make_passband(filter_name=filter_name)

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
            filter_curve = make_passband(filter_file=filter_file)
        else:
            filter_curve = make_passband(filter_name=filter_name)

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
            filter_curve = make_passband(filter_file=filter_file)
        else:
            filter_curve = make_passband(filter_name=filter_name)

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
    Copied from SimCADO
    Retrieve the Vega spectrum from stsci and return it in synphot format

    Notes
    -----
    To access wavelength and fluxes use::

        wave, flux = vega_sp._get_arrays(wavelengths=None)

    """
    location = "http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_009.fits"
    remote = synphot.specio.read_remote_spec(location, cache=True)
    header = remote[0]
    wave = remote[1]
    flux = remote[2]
    url = 'Vega from ' + location
    meta = {'header': header, 'expr': url}
    vega_sp = Spextrum(modelclass=Empirical1D, points=wave, lookup_table=flux, meta=meta)
    return vega_sp




