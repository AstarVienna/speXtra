# -*- coding: utf-8 -*-
"""
Spyctra: A python tool to manage and manipulate astronomical spectra

returned spectra should  be in synphot format: do not alter attributes

Database is implemented as a yaml file containing the basic information of the libraries
each library is described by a yaml file

Classes for describing the database and libraries are being implemented
TODO: Consider @dataclass for better and more concise description of these classes (only python 3.7 though)
"""
import numbers
import warnings

import numpy as np

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.utils.data import download_file
from astropy.constants import c

import synphot
from synphot import (units, SourceSpectrum, SpectralElement, Observation, BaseUnitlessSpectrum)
from synphot.models import (Empirical1D, GaussianFlux1D, Box1D, ConstFlux1D)

import tynt

from .database import get_template

# default filter definitions now in data/default_filters.yml, including GALEX, Spitzer, HST and Y-band filters
# with open("data/default_filters.yml") as f:
#    FILTER_DEFAULTS = yaml.safe_load(f)


class Spectrum(SourceSpectrum):
    """
    Class to handle spectra. This class stores and manipulates the spectra.

    This class can be initialized with a remote file which will be downloaded from
    the database or with a synphot.Spectrum

    Parameters
    ----------
    template_name: Name of the template to download with format library/template
    modelclass, kwargs
        See `BaseSpectrum`.

    """

    def __init__(self, template_name=None, modelclass=None, **kwargs):

        self.template_name = template_name

        if self.template_name is not None:
            meta, lam, flux = self.__loader()
            modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)
        if modelclass is not None:
            modelclass = modelclass
        else:
            raise ValueError("please define a spectra")

        super().__init__(modelclass, **kwargs)

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
        self.resolution = meta["resolution"]
        self.wave_unit = meta["wave_unit"]
        self.flux_unit = meta["flux_unit"]
        self.wave_column_name = meta["wave_column_name"]
        self.flux_column_name = meta["flux_column_name"]
        self.ile_extension = meta["file_extension"]
        # This should be atrributes too

        # make try and except here to catch most problems
        if data_type == "fits":
            meta, lam, flux = synphot.specio.read_fits_spec(location)
        else:
            meta, lam, flux = synphot.specio.read_ascii_spec(location)

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
        filename
        format: format of the spectra accepted by specutils (see specutils documentation)

        Returns
        -------


        """
        try:
            from specutils import Spectrum1D
        except ImportError as ie:
            print(ie, "specutils not installed, cannot import spectra")

        spec1d = Spectrum1D.read(filename, format=format, **kwargs)
        meta = {'header': dict(fits.getheader(filename))}
        lam = spec1d.spectral_axis.value * u.AA
        flux = spec1d.flux.value * units.FLAM

        modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)
        return Spectrum(modelclass=modelclass)

    def redshift(self, z=0, vel=0):
        """
        Redshift or blueshift a spectra

        Parameters
        ----------
        spectrum: a synphot spectra
        z: redshift
        vel: radial velocity,  if no unit are present it is assumed to be in m/s

        Returns
        -------
        a synphot SourceSpectrum

        """
        if vel != 0:
            if isinstance(vel, u.Quantity) is False:
                vel = vel * u.m / u.s # assumed to be in m/s

            z = (vel.to(u.m / u.s) / c).value
        if z <= -1:
            raise ValueError("Redshift or velocity unphysical")

        lam = self.model.points[0] * (1 + z)
        flux = self.model.lookup_table
        meta = self.model.meta

        modelclass = SourceSpectrum(Empirical1D, points=lam, lookup_table=flux, meta=meta)
        return Spectrum(modelclass=modelclass)

    @classmethod
    def rebin_spectra(cls, new_waves):
        """
        Rebin a synphot spectra to a new wavelength grid conserving flux.
        Grid does not need to be linear and can be at higher or lower resolution

        TODO: To resample the spectra at lower resolution a convolution is first needed. Implement!
        TODO: Return the new spectra in the input wavelengths units
        TODO: Check this!!

        Parameters
        ----------
        spectra: a synphot spectra
        new_waves: an array of the output wavelenghts in Angstroms but other units can be
            specified

        Returns
        -------

        A synphot spectra in the new wavelengths

        """
        if isinstance(new_waves, u.Quantity):
            new_waves = new_waves.to(u.AA).value

        waves = cls.waveset.value
        f = np.ones(len(waves))
        filt = SpectralElement(Empirical1D, points=waves, lookup_table=f)
        obs = Observation(cls, filt, binset=new_waves, force='taper')
        newflux = obs.binflux
        rebin_spec = SourceSpectrum(Empirical1D, points=new_waves, lookup_table=newflux, meta=cls.meta)

        return Spectrum(modelclass=rebin_spec)

    def add_emission_line(self, center, flux, fwhm):
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

        for c, x, f in zip(center, flux, fwhm):
            g_em = SourceSpectrum(GaussianFlux1D(mean=c, total_flux=x, fwhm=f))
            sp = self.__class__(self.model + g_em.model)

        return sp

    def add_absorption_line(self, center, ew, fwhm):
        """
        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)
        Add a absorption line of to a spectrum with center, fwhm and equivalent width specified by the user
        It also supports emission lines if ew is negative

        Parameters
        ----------
        spectrum
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

        for c, e, f in zip(center, ew, fwhm):
            sign = -1 * np.sign(e)  # to keep the convention that EL are negative and ABS are positive
            left, right = center - np.abs(e / 2), center + np.abs(e / 2)
            wavelengths = self.waveset[(self.waveset.value >= left) & (self.waveset.value <= right)]
            fluxes = units.convert_flux(wavelengths=wavelengths, fluxes=self(wavelengths),
                                        out_flux_unit=units.FLAM)
            flux = np.trapz(fluxes.value, wavelengths.value)
            g_abs = Spectrum(GaussianFlux1D(total_flux=sign * flux, mean=c, fwhm=f))
            sp = self.__class__(self.model + g_abs.model)
            if (sp(wavelengths).value < 0).any():
                warnings.warn("Warning: Flux<0 for specified EW and FHWM, setting it to Zero")
                waves = sp.waveset[sp(sp.waveset) < 0]
                zero_sp = SourceSpectrum(Empirical1D, points=waves, lookup_table=-1 * sp(waves).value)
                sp = self.__class__(sp.model + zero_sp.model)

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
        sp = self.__class__(self.model * extinction.model)

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
        sp = self.__class__(self.model / extinction.model)

        return sp

    def smooth(self, kernel):
        pass

    @classmethod
    def zero_mag_spectrum(cls, system_name="AB"):
        """
        Creates a spectrum with zero magnitude in the preferred system
        Parameters
        ----------
        system_name: AB, Vega or ST, default AB

        Returns
        -------

        """
        mag = 0
        if system_name.lower() in ["vega"]:
            vega = get_vega_spectrum()
            spec = vega * 10 ** (-0.4 * mag)  # is this necessary?
        elif system_name.lower() in ["ab"]:
            spec = SourceSpectrum(ConstFlux1D, amplitude=mag * u.ABmag)
        elif system_name.lower() in ["st", "hst"]:
            spec = SourceSpectrum(ConstFlux1D, amplitude=mag * u.STmag)

        return cls(modelclass=spec)

    @classmethod
    def black_body_spectrum(cls, temperature, wmin, wmax):
        """
        Produce a blackbody spectrum for a given temperature and cut it between wmin and wmax
        Parameters
        ----------
        temperature
        wmin
        wmax

        Returns
        -------

        """
        pass

    def scale_to_magnitude(self, magnitude, unit):
        pass

    def _validate_other_add_sub(self, other):
        """
        Conditions for other to satisfy before add/sub.
        Copied from SourceSpectrum so additions can happen with SourceSpectrum
        """
        if not isinstance(other, (self.__class__, SourceSpectrum)):
            raise exceptions.IncompatibleSources(
                'Can only operate on {0}.'.format(self.__class__.__name__))

    def _validate_other_mul_div(self, other):
        """Conditions for other to satisfy before mul/div."""
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


#------------------ BEGIN ---------------------------------------------
# This has been copied from scopesim.effects.ter_curves_utils.py

def download_svo_filter(filter_name):
    """
    Query the SVO service for the true transmittance for a given filter
    Copied 1 to 1 from tynt by Brett Morris
    Parameters
    ----------
    filter_name : str
        Name of the filter as available on the spanish VO filter service
        e.g: ``Paranal/HAWKI.Ks``
    Returns
    -------
    filt_curve : ``synphot.SpectralElement``
        Astronomical filter object.
    """
    path = download_file('http://svo2.cab.inta-csic.es/'
                         'theory/fps3/fps.php?ID={}'.format(filter_name),
                         cache=True)

    true_transmittance = Table.read(path, format='votable')
    wave = true_transmittance['Wavelength'].data.data * u.Angstrom
    trans = true_transmittance['Transmission'].data.data
    filt_curve = SpectralElement(Empirical1D, points=wave, lookup_table=trans)

    return filt_curve


def get_filter(filter_name):
    # first check locally
    # check generics
    # check spanish_vo
    path = find_file(filter_name, silent=True)

    if path is not None:
        tbl = ioascii.read(path)
        wave = quantity_from_table("wavelength", tbl, u.um).to(u.um)
        filt = SpectralElement(Empirical1D, points=wave,
                               lookup_table=tbl["transmission"])
    elif filter_name in FILTER_DEFAULTS:
        filt = download_svo_filter(FILTER_DEFAULTS[filter_name])
    else:
        try:
            filt = download_svo_filter(filter_name)
        except:
            filt = None

    return filt





def zero_mag_flux(filter_name, photometric_system, return_filter=False):
    """
    Returns the zero magnitude photon flux for a filter
    Acceptable filter names are those given in
    ``scopesim.effects.ter_curves_utils.FILTER_DEFAULTS`` or a string with an
    appropriate name for a filter in the Spanish-VO filter-service. Such strings
    must use the naming convention: observatory/instrument.filter. E.g:
    ``paranal/HAWKI.Ks``, or ``Gemini/GMOS-N.CaT``.
    Parameters
    ----------
    filter_name : str
        Name of the filter - see above
    photometric_system : str
        ["vega", "AB", "ST"] Name of the photometric system
    return_filter : bool, optional
        If True, also returns the filter curve object
    Returns
    -------
    flux : float
        [PHOTLAM]
    filt : ``synphot.SpectralElement``
        If ``return_filter`` is True
    """

    filt = get_filter(filter_name)
    spec = get_zero_mag_spectrum(photometric_system)

    obs = Observation(spec, filt)
    flux = obs.effstim(flux_unit=PHOTLAM)

    if return_filter:
        return flux, filt
    else:
        return flux


def scale_spectrum(spectrum, filter_name, amplitude):
    """
    Scales a SourceSpectrum to a value in a filter
    Parameters
    ----------
    spectrum : synphot.SourceSpectrum
    filter_name : str
        Name of a filter from
        - a local instrument package (available in ``rc.__search_path__``)
        - a generic filter name (see ``ter_curves_utils.FILTER_DEFAULTS``)
        - a spanish-vo filter service reference (e.g. ``"Paranal/HAWKI.Ks"``)
    amplitude : ``astropy.Quantity``, float
        The value that the spectrum should have in the given filter. Acceptable
        astropy quantities are:
        - u.mag : Vega magnitudes
        - u.ABmag : AB magnitudes
        - u.STmag : HST magnitudes
        - u.Jy : Jansky per filter bandpass
        Additionally the ``FLAM`` and ``FNU`` units from ``synphot.units`` can
        be used when passing the quantity for ``amplitude``:
    Returns
    -------
    spectrum : synphot.SourceSpectrum
        Input spectrum scaled to the given amplitude in the given filter
    Examples
    --------
    ::
        >>> from scopesim.effects.ter_curves_utils as ter_utils
        >>>
        >>> spec = ter_utils.vega_spectrum()
        >>> vega_185 = ter_utils.scale_spectrum(spec, "Ks", -1.85 * u.mag)
        >>> ab_0 = ter_utils.scale_spectrum(spec, "Ks", 0 * u.ABmag)
        >>> jy_3630 = ter_utils.scale_spectrum(spec, "Ks", 3630 * u.Jy)
    """

    if isinstance(amplitude, u.Quantity):
        if amplitude.unit.physical_type == "spectral flux density":
            if amplitude.unit != u.ABmag:
                amplitude = amplitude.to(u.ABmag)
            ref_spec = ab_spectrum(amplitude.value)

        elif amplitude.unit.physical_type == "spectral flux density wav":
            if amplitude.unit != u.STmag:
                amplitude = amplitude.to(u.STmag)
            ref_spec = st_spectrum(amplitude.value)

        elif amplitude.unit == u.mag:
            ref_spec = vega_spectrum(amplitude.value)

        else:
            raise ValueError("Units of amplitude must be one of "
                             "[u.mag, u.ABmag, u.STmag]: {}".format(amplitude))
    else:
        ref_spec = vega_spectrum(amplitude)

    filt = get_filter(filter_name)
    ref_flux = Observation(ref_spec, filt).effstim(flux_unit=PHOTLAM)

    real_flux = Observation(spectrum, filt).effstim(flux_unit=PHOTLAM)
    scale_factor = ref_flux / real_flux
    spectrum *= scale_factor.value

    return spectrum




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
    location = "http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_008.fits"
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


def get_filter(name=None, filename=None, wave_unit=u.AA):
    """
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

    """

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


def photons_in_range(spectra, wave_min, wave_max, area, bandpass=None):
    """
    Return the number of photons between wave_min and wave_max or within
    a bandpass (filter)

    TODO: Write wrapper functions make_synphot_bandpass and make_synphot_spectra
        to allow a variety of bandpasses and spectra.



    Parameters
    ----------
    spectra: a synphot spectrum
    wave_min
        [Angstrom]
    wave_max
        [Angstrom]
    area : Quantity
        [cm2]
    bandpass : SpectralElement


    Returns
    -------
    counts : u.Quantity array

    """
    if isinstance(area, u.Quantity):
        area = area.to(u.cm**2).value  #
    if isinstance(wave_min, u.Quantity):
        wave_min = wave_min.to(u.Angstrom).value
    if isinstance(wave_max, u.Quantity):
        wave_max = wave_max.to(u.Angstrom).value
    if isinstance(spectra, list) is False:
        spectra = [spectra]

    if bandpass is None:
        # this makes a bandpass out of wmin and wmax
        mid_point = 0.5*(wave_min + wave_max)
        width = abs(wave_max - wave_min)
        bandpass = SpectralElement(Box1D, amplitude=1, x_0=mid_point, width=width)

    if (bandpass is not None) and (isinstance(bandpass, synphot.spectrum.SpectralElement) is False) :
        # bandpass = make_synphot_bandpass(bandpass) # try to make a synphot bandpass from e.g. filter file
        pass

    counts = []
    for spec in spectra:
        if isinstance(spec, synphot.spectrum.SourceSpectrum) is False:
            #spec = make_synphot_spectrum(spec) # Try to make a synphot spectrum from e.g. file/np.array
            pass

        obs = Observation(spec, bandpass)
        counts.append(obs.countrate(area=area*u.m**2).value)

    counts = np.array(counts) * u.ph * u.s**-1

    return counts



def black_body(t, wmin, wmax):
    """
    TODO: It is needed?
    Unitility function to create SourceSpectrum with a black body with temperature T

    Parameters
    ----------
    t: Temperature in Kelvin
    wmin
    wmax

    Returns
    -------


    """
    pass




