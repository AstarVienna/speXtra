# -*- coding: utf-8 -*-
"""
Spyctra: A python tool to manage and manipulate astronomical spectra

returned spectra should  be in synphot format: do not altere attributes


Probably we need functions to parse particular spectra, like Pickles, Brown, Kinney-Calzetti, etc
Another possibility is to download all these catalogs and put them in a single format.


We need to create a database of remote spectra like this:

collection | name | type   | flux_unit | wave_unit | location | wave_col | flux_col

Kinney     | Sba1 | galaxy | FLAM      | Angstrom  |  URL   |
Pickles    | A0V  | AGN    | etc       |  um       |     |
etc

adding spectra should be as easy as table.add_row()

What about also a user defined database of spectra?

two options

1.- Add the spectra to the database allowing also paths besides URLs, and collection='User'

2.- A different database containing the same info

3.- A different database containing arrays that can be read as wavelengths and fluxes
    name | type | lambda      | flux         | flux_unit | wave_unit
    sp1  | gal1 | (3000,4000) | (1e-3, 1e-4)
    sp2  | gal2


"""

import numpy as np
import synphot
from synphot import units
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.utils.data import download_file
from astropy.constants import c
import tynt
import warnings


class Spectra:
    """
    This class stores and manipulates the spectra. Whenever applies, it returns a synphot spectra

    TODO: Write the methods first as individual functions and incorporate later
    """

    def from_file(self, filename):
        pass

    def from_database(self, name):
        pass

    def plot(self, wave_units=u.angstrom, flux_unit="flam"):

        wavelenghts = self.waveset.to(wave_unit, equivalencies=u.spectral())
        self.plot(wavelenghts=wavelenghts, flux_unit=flux_unit)
        pass

    def redshift(self, z=0, vel=0):
        pass

    def attenuate(self, curve):
        pass

    def rebin(self, array):
        pass

    def smooth(self, kernel):
        pass

    def scale_to_magnitude(self, magnitude, unit):
        pass



def load_spectra(filename, wave_unit=u.AA, flux_unit=units.FLAM, wave_col='WAVELENGTH', flux_col='FLUX', ext=1):
    """
    This function try to load a spectra and return it as a synphot.SourceSpectrum object.
    It will attempt to load the spectra from different format but it might nevertheless fail.

    Parameters
    ----------
    filename: A file containing the 1D spectra
    wave_unit: optional,
                wavelength units, either synphot.units or astropy.units
    flux_unit: optional,
                flux units, either synphot.units or astropy.units
    wave_col: optional,
                The column name containing the wavelengths if filename is a fits table, default WAVELENGTH
    flux_col: optional,
               The column name containing the wavelengths if filename is a fits table, default FLUX
    ext: optional,
         The extension number where the spectra is located if file is a fits table or image, default 1


    Returns
    -------

    A synphot.SourceSpectrum

    """
    # This should load most of the tables
    try:
        if filename.lower().endswith("fit") or filename.lower().endswith("fits"):
            sp = synphot.SourceSpectrum.from_file(filename=filename, wave_unit=wave_unit, flux_unit=flux_unit,
                                          wave_col=wave_col, flux_col=flux_col, ext=ext)
        else:
            sp = synphot.SourceSpectrum.from_file(filename=filename, wave_unit=wave_unit, flux_unit=flux_unit)

    except AttributeError as e:
     # Try to load the 1D fits images-spectra here
        print(e)
     #   try:
     #       sp = load_1dfits(filename,                            wave_unit, flux_unit, wave_col, flux_col, ext)

     #   except AttributeError as e:
     #       print(e)

    return sp



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
    vega_sp = synphot.SourceSpectrum(synphot.models.Empirical1D,
                                     points=wave, lookup_table=flux, meta=meta)
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

    bandpass = synphot.SpectralElement(synphot.models.Empirical1D(points=waves,
                                                                  lookup_table=trans))

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
        [um]
    wave_max
        [um]
    area : Quantity
        [m2]
    bandpass : SpectralElement


    Returns
    -------
    counts : u.Quantity array

    """
    if isinstance(area, u.Quantity):
        area = area.to(u.m**2).value  # if not unit is given, area is assumed in m2

    if isinstance(wave_min, u.Quantity):
        wave_min = wave_min.to(u.Angstrom).value
    else:
        wave_min *= 1E4

    if isinstance(wave_max, u.Quantity):
        wave_max = wave_max.to(u.Angstrom).value
    else:
        wave_max *= 1E4  # if not unit is given, wavelength is assumed in Angstrom

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

        obs = synphot.Observation(spec, bandpass)
        counts.append(obs.countrate(area=area*u.m**2).value)

    counts = np.array(counts) * u.ph * u.s**-1

    return counts



def rebin_spectra(spectra, new_waves):
    """
    Rebin a synphot spectra to a new wavelength grid conserving flux.
    Grid does not need to be linear and can be at higher or lower resolution

    TODO: To resample the spectra at lower resolution a convolution is first needed. Implement!
    TODO: Return the new spectra in the input wavelengths

    Parameters
    ----------
    spectra: a synphot spectra
    new_waves: an array of the output wavelenghts in Angstroms but other units can be
        specified


    Returns
    -------

    A synphot spectra in the new wavelengths

    """
    if isinstance(spectra, synphot.spectrum.SourceSpectrum) is False:
        # spec = make_synphot_spectrum(spec) # Try to make a synphot spectrum from e.g. file/np.array
        raise ValueError(spec, "is not a synphot spectra!")

    if spectra.waveset is None:
        raise ValueError("spectra doesn't have a defined waveset")

    if isinstance(new_waves, u.Quantity):
        new_waves = new_waves.to(u.Angstrom).value

    waves = spectra.waveset.value
    f = np.ones(len(waves))

    filt = synphot.SpectralElement(synphot.models.Empirical1D,
                                   points=waves, lookup_table=f)
    obs = synphot.Observation(spectra, filt, binset=new_waves, force='taper')

    newflux = obs.binflux

    rebin_spec = synphot.SourceSpectrum(synphot.models.Empirical1D,
                                        points=new_waves, lookup_table=newflux, meta=spectra.meta)

    return rebin_spec


def scale_to_magnitude(spectra, magnitude, passband, units=u.ABmag):
    """
    Scale spectra to a magnitude

    Parameters
    ----------
    spectra
    magnitude
    units

    Returns
    -------
    a scale synphot spectra
    """
    pass


def get_magnitude(spectra, passband, units=u.ABmag):
    """
    get magnitude in the passband from spectra spectra to a magnitude

    Parameters
    ----------
    spectra
    magnitude
    units

    Returns
    -------
    a magnitude in the specified units
    """
    pass


def redshift(spectrum, z=0, vel=0):
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
    if z != 0:
        spectrum.z = z

    if vel != 0:

        if isinstance(vel, u.Quantity) is False:

            vel = vel * u.m / u.s

        z = vel.to(u.m/u.s)/c
        spectrum.z = z.value

    return spectrum


def attenuate(spectrum, curve, AV=None, EBV=0):
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

    extinction = synphot.ReddeningLaw.from_extinction_model(curve).extinction_curve(EBV)
    sp = spectrum * extinction

    return sp


def deredden(spectrum, curve, EBV=0):
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

    extinction = synphot.ReddeningLaw.from_extinction_model(curve).extinction_curve(EBV)
    sp = spectrum / extinction

    return sp







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


def add_emission_line(spectrum, center, flux, fwhm):
    """

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
        flux = flux.to(u.erg/(u.cm**2 * u.s)).value
    if isinstance(fwhm, u.Quantity) is True:
        fwhm = fwhm.to(u.AA).value

    g_em = synphot.SourceSpectrum(synphot.models.GaussianFlux1D(total_flux=flux,
                                                                mean=center, fwhm=fwhm))
    sp = spectrum + g_em

    return sp


def add_absorption_line(spectrum, center, ew, fwhm):
    """
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

    sign = -1 * np.sign(ew)  # to keep the convention that EL are negative and ABS are positive
    left, right = center - np.abs(ew/2), center + np.abs(ew/2)
    wavelengths = spectrum.waveset[(spectrum.waveset.value >= left) & (spectrum.waveset.value <= right)]

    fluxes = synphot.units.convert_flux(wavelengths=wavelengths, fluxes=spectrum(wavelengths),
                                        out_flux_unit=units.FLAM)
    flux = np.trapz(fluxes.value, wavelengths.value)

    g_abs = synphot.SourceSpectrum(synphot.models.GaussianFlux1D(total_flux=sign * flux,
                                                                 mean=center, fwhm=fwhm))
    sp = spectrum + g_abs

    if (sp(wavelengths).value < 0).any():
        warnings.warn("Warning: Flux<0 for specified EW and FHWM, setting it to Zero")
        waves = sp.waveset[sp(sp.waveset) < 0]
        zero_sp = synphot.SourceSpectrum(synphot.models.Empirical1D,
                                         points=waves, lookup_table=-1 * sp(waves).value)
        sp = sp + zero_sp

    return sp





def list_spectra(keyword):
    """
    List the spectra available for download
    """
    pass


def get_spectra(keyword):
    """
    This should download the spectra from the interwebs
    """
    pass


def get_pickles(type):
    pass




def load_1dfits(filename, wave_unit=u.AA, flux_unit=units.FLAM, wave_col='WAVELENGTH', flux_col='FLUX', ext=1):
    """
    This function load spectra stored as fits images and return them in fits format

    TODO: MUCH TO DO HERE.

    Parameters
    ----------
    filename: A file containing the 1D spectra
    wave_unit: optional,
                wavelength units, either synphot.units or astropy.units
    flux_unit: optional,
                flux units, either synphot.units or astropy.units
    wave_col: optional,
                The column name containing the wavelengths if filename is a table
    flux_col: optional,
               The column name containing the wavelengths if filename is a table
    ext: optional,
         The extension number where the spectra is located if file is a fits table or image, default 1


    Returns
    -------

    A synphot.SourceSpectrum

    """

    hdu = fits.open(filename)

    header_primary = hdu[0].header
    data = hdu[ext].data
    if len(hdu) > 1:
        header_ext = hdu[ext].header
        data = hdu[ext].data
    else:
        header_ext = header_primary

    hdu.close()
    try:
        crval1 = header_ext["CRVAL1"]
        crpix1 = header_ext["CRPIX1"]
        cdelt1 = header_ext["CDELT1"]
        #cunit1 = header_ext["CUNIT1"]
        bunit  = header_ext["BUNIT"]

    except KeyError as e:
        print(e, "keys not found")
        try:
            crval1 = header_primary["CRVAL1"]
            crpix1 = header_primary["CRPIX1"]
            cdelt1 = header_primary["CDELT1"]
         #   cunit1 = header_primary["CUNIT1"]
            bunit = header_primary["BUNIT"]

        except KeyError as e:
            print(e, "keys not found")

    try:
        spl_units = bunit.split(" ")
        factor = eval(spl_units.pop(0))
        unit = u.Unit(" ".join(spl_units))

    except ValueError as e:
        print(e, "Units not recognized")

    wavelengths = crval1 + np.arange(crpix1 - 1, len(data), 1)*cdelt1
    values = data
    fluxes = values*factor*unit
    sp = synphot.SourceSpectrum(synphot.models.Empirical1D, points=wavelengths, lookup_table=fluxes,
                        meta={'header': header_primary})

    return sp


