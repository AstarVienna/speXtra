# -*- coding: utf-8 -*-
"""
speXtra: A python tool to manage and manipulate astronomical spectra
"""

import numbers
import warnings

import numpy as np

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.constants import c
from astropy.modeling.models import Scale

import synphot
from synphot import (units, SourceSpectrum, SpectralElement, Observation, BaseUnitlessSpectrum)
from synphot.models import (Empirical1D, GaussianFlux1D, Box1D, ConstFlux1D, BlackBody1D)
from synphot import exceptions

import tynt

from .database import get_template, get_filter


def make_passband(filter_name=None, filter_file=None, wave_unit=u.Angstrom):
    """
    Make a SpectralElement (synphot passband) from user specified filter in the database.
    Optionally, specify a file in disk and the wavelength units.

    Parameters
    ----------
    filter_name: str, a filter name expressed as ``instrument/filter_name``
    filter_file: Optionally, make a pasband from a local file
    wave_unit: astropy.unit

    Returns
    -------
    passband: a synphot.SpectralElement
    """
    if filter_file is not None:
        try:
            meta, wave, trans = synphot.specio.read_spec(filter_file,
                                                         wave_unit=wave_unit,
                                                         flux_unit='transmission')
        except (FileNotFoundError, synphot.exceptions.SynphotError) as e:
            print("File not found or malformed", e)
    else:
        path, meta = get_filter(filter_name)
        if meta is None:  # it's a svo filter
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
    the database or with a synphot.Spectrum


    Parameters
    ----------
    template_name: Name of the template to download with format library/template e.g. "kc96/s0
    modelclass, kwargs
        See `BaseSpectrum`.

    """

    def __init__(self, template_name=None, modelclass=None, **kwargs):

        self.template_name = template_name

        if self.template_name is not None:
            meta, lam, flux = self.__loader
            modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)
        if modelclass is not None:
            modelclass = modelclass
        else:
            raise ValueError("please define a spectra")

        super().__init__(modelclass, **kwargs)

    @property
    def __loader(self):
        """
        Load a template from the database

        TODO: Checks for units etc in the library to correctly call read_*__spec for most possible cases
        Parameters
        ----------
        template_name: The name of the spectral template in the speclibrary, format: library/template_name

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux

        """
        location, meta = get_template(self.template_name)
        print(location)

        data_type = meta["data_type"]
        try:  # it should also try to read it from the file directly
            self.wave_unit = units.validate_unit(meta["wave_unit"])
            self.flux_unit = units.validate_unit(meta["flux_unit"])
        except synphot.exceptions.SynphotError:
            self.wave_unit = u.AA
            self.flux_unit = units.FLAM

        self.resolution = meta["resolution"] * self.wave_unit
        self.wave_column_name = meta["wave_column_name"]  # same here
        self.flux_column_name = meta["flux_column_name"]
        self.file_extension = meta["file_extension"]

        # make try and except here to catch most problems
        if data_type == "fits":
            meta, lam, flux = synphot.specio.read_fits_spec(location, ext=1,
                                                            wave_unit=self.wave_unit,
                                                            flux_unit=self.flux_unit,
                                                            wave_col=self.wave_column_name,
                                                            flux_col=self.flux_column_name)
        else:
            meta, lam, flux = synphot.specio.read_ascii_spec(location,
                                                             wave_unit=self.wave_unit, flux_unit=self.flux_unit)

        self.wmin = np.min(lam)
        self.wmax = np.max(lam)

        return meta, lam, flux

    @classmethod
    def from1dspec(cls, filename, format="wcs1d-fits", **kwargs):
        """
        This function _tries_ to create a Spectrum from 1d fits files.
        It relies in specutils for its job

        TODO: More checks from units, etc.
        Parameters
        ----------
        filename: str a filename with the spectra
        format: format of the spectra accepted by specutils (see specutils documentation)

        Returns
        -------
        a Spextrum instance

        """
        try:
            from specutils import Spectrum1D
        except ImportError as ie:
            print(ie, "specutils not installed, cannot import spectra")

        spec1d = Spectrum1D.read(filename, format=format, **kwargs)
        meta = {'header': dict(fits.getheader(filename))}
        lam = spec1d.spectral_axis.value * u.AA
        flux = spec1d.flux.value * units.FLAM

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
                vel = vel * u.m / u.s # assumed to be in m/s

            z = (vel.to(u.m / u.s) / c).value
        if z <= -1:
            raise ValueError("Redshift or velocity unphysical")

        lam = self.waveset * (1 + z)
        flux = self(self.waveset)
        meta = self.meta
        modelclass = SourceSpectrum(Empirical1D,
                                    points=lam, lookup_table=flux, meta=meta)

        return Spextrum(modelclass=modelclass)

    def rebin_spectra(self, new_waves):
        """
        Rebin a synphot spectra to a new wavelength grid conserving flux.
        Grid does not need to be linear and can be at higher or lower resolution

        TODO: To resample the spectra at lower resolution a convolution is first needed. Implement!

        Parameters
        ----------
        new_waves: an array of the output wavelenghts in Angstroms but other units can be
            specified

        Returns
        -------

        A synphot spectra in the new wavelengths

        """
        if isinstance(new_waves, u.Quantity):
            new_waves = new_waves.to(u.AA).value

        waves = self.waveset.value   # else assumed to be angstroms
        f = np.ones(len(waves))
        filt = SpectralElement(Empirical1D, points=waves, lookup_table=f)
        obs = Observation(self, filt, binset=new_waves, force='taper')
        newflux = obs.binflux
        rebin_spec = SourceSpectrum(Empirical1D, points=new_waves, lookup_table=newflux, meta=self.meta)

        return Spextrum(modelclass=rebin_spec)

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
        the spectrum with the emission line

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
        #print(center, flux, fwhm)

        sp = self.__class__(modelclass=self.model)
        for c, x, f in zip(center, flux, fwhm):
            g_em = SourceSpectrum(GaussianFlux1D(mean=c, total_flux=x, fwhm=f))

            sp = self.__class__(modelclass=sp.model + g_em.model)

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

        sp = self.__class__(modelclass=self.model)
        for c, e, f in zip(center, ew, fwhm):
            sign = -1 * np.sign(e)  # to keep the convention that EL are negative and ABS are positive
            left, right = center - np.abs(e / 2), center + np.abs(e / 2)
            wavelengths = self.waveset[(self.waveset.value >= left) & (self.waveset.value <= right)]
            fluxes = units.convert_flux(wavelengths=wavelengths,
                                        fluxes=self(wavelengths),
                                        out_flux_unit=units.FLAM)
            flux = np.trapz(fluxes.value, wavelengths.value)
            g_abs = SourceSpectrum(GaussianFlux1D(total_flux=sign * flux, mean=c, fwhm=f))
            sp = self.__class__(modelclass=sp.model + g_abs.model)

            if (sp(wavelengths).value < 0).any():
                warnings.warn("Warning: Flux<0 for specified EW and FHWM, setting it to Zero")
                waves = sp.waveset[sp(sp.waveset) < 0]
                zero_sp = SourceSpectrum(Empirical1D, points=waves, lookup_table=-1 * sp(waves).value)
                sp = self.__class__(modelclass=sp.model + zero_sp.model)

        return sp

    def redden(self, curve, Av=None, Ebv=0, Rv=3.1):
        """
        This function attenuate  a spectrum with a extinction curve normalized to a E(B-V)

        TODO: We also need a library of extinction curves and a way to read them

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
            Ebv = Rv * Ebv

        extinction = synphot.ReddeningLaw.from_extinction_model(curve).extinction_curve(Ebv)
        sp = self.__class__(modelclass=self.model * extinction.model)

        return sp

    def deredden(self, curve, Av=None, Ebv=0, Rv=3.1):
        """
        This function de-redden a spectrum.

        TODO: We also need a library of extinction curves and a way to read them

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
            Ebv = Rv * Ebv

        extinction = synphot.ReddeningLaw.from_extinction_model(curve).extinction_curve(Ebv)
        sp = self.__class__(modelclass=self.model / extinction.model)

        return sp

    def smooth(self, kernel):
        pass

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
        if wavelengths is None: # set a default waveset
            wavelengths = np.geomspace(100, 5e4, num=5000) * u.AA # constant R~800

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
        a black-body spectrum
        """
        bb = SourceSpectrum(BlackBody1D, temperature=temperature)
        sp = cls(modelclass=bb)

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
        real_flux = Observation(SourceSpectrum(modelclass=self),
                                filter_curve).effstim(flux_unit=units.PHOTLAM)
        mag = -2.5*np.log10(real_flux.value/ref_flux.value)

        return mag * unit

    def photons_in_range(self, wmin, wmax, area=1*u.cm**2,
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

    def add_noise(self):
        """
        Create a spectra with a given S/N in the wavelength range
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
    vega_sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, meta=meta)
    return vega_sp


def get_filter_systems():
    """
    Return a set of the different filter system available

    Returns
    -------

    """
    filters = tynt.FilterGenerator().available_filters()
    systems = {f.split("/")[0] for f in filters}
    return systems


def get_filter_names(system=None):
    """
    This function just returns the filters available from tynt
    if system= None returns all

    Returns
    -------

    """
    filter_list = tynt.FilterGenerator().available_filters()
    ord_list = [[f for f in filter_list if s in f] for s in get_filter_systems()]
    flat_list = [item for sublist in ord_list for item in sublist]

    if system is not None:
        flat_list = [f for f in filter_list if s in f]

    return flat_list


"""
def get_filter(name=None, filename=None, wave_unit=u.AA):

    Return a synphot SpectralElement (bandpass) from a filter in the SVO Filter Profile Service
    or from a local file (only ascii is supported atm)

    TODO: It should be also able to read a SpectralElement (bandpass) from synphot
    See ScopeSim.effects.TER_curve_utils


    Parameters
    ----------
    name: Name of the filter in the SVO Filter Profile Service, an incomplete list supplied by tynt
          can be obtained with get_filter_names()
    filename: Filename where filter is stored, only ascii is supported, col1 is assumed to be wavelength, col2 transmittance
    wave_unit: unit of the wavelength column, default is u.AA (Angstroms)

    Returns
    -------



    if name is not None:
        try:
            f = tynt.FilterGenerator()
            filt = f.download_true_transmittance(name)
            waves = filt.wavelength.value
            trans = filt.transmittance

        except ValueError as e:
            print(e, "Filter not found in SVO Filter Profile Service")

    if filename is not None:
        try:
            filter_table = Table.read(filename, format="ascii")
            waves = filter_table["col1"].data * wave_unit
            waves = waves.to(u.AA).value
            trans = filter_table["col2"].data

        except FileNotFoundError as e:
            print(e, "File not found")

    bandpass = SpectralElement(Empirical1D(points=waves, lookup_table=trans))

    return bandpass

"""



