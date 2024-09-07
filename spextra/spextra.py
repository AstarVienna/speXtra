# -*- coding: utf-8 -*-
"""speXtra: A python tool to manage and manipulate astronomical spectra."""

import numbers
import warnings
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from more_itertools import always_iterable

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.modeling.models import Scale

from synphot import (units, SourceSpectrum, SpectralElement, Observation,
                     BaseUnitlessSpectrum, ReddeningLaw, utils)
from synphot.models import (Empirical1D, GaussianFlux1D, Gaussian1D, Box1D,
                            ConstFlux1D, BlackBody1D, PowerLawFlux1D)
from synphot.specio import read_ascii_spec, read_fits_spec
from synphot import exceptions

from .database import DEFAULT_DATA
from .containers import SpectrumContainer, FilterContainer, ExtCurveContainer
from .downloads import download_svo_filter
from .exceptions import SpextraError, ArgumentError, ConstructorError


__all__ = ["Spextrum", "Passband", "ExtinctionCurve", "get_vega_spectrum"]


class Passband(SpectralElement, FilterContainer):
    """
    Class to handle astronomical filters.

    TODO: Implement proper __add__, __sub__, __mul__, etc

    """

    DEFAULT_FILTERS = DEFAULT_DATA.filters

    def __init__(self, filter_name=None, modelclass=None, **kwargs):
        if filter_name is not None:
            self.from_filter_name(filter_name)
        elif modelclass is not None:
            self.from_modelclass(modelclass, **kwargs)
        else:
            raise ConstructorError(
                "Either filter_name or modelclass must be passed."
            )

    # @classmethod
    def from_filter_name(self, filter_name: str):
        filter_name = self.DEFAULT_FILTERS.get(filter_name, filter_name)

        FilterContainer.__init__(self, filter_name)

        if self.is_in_library and self.library.is_in_database:
            self._database_loader()
        else:
            try:  # try to download it from SVO
                self._svo_loader(filter_name)
            except ValueError:
                raise SpextraError(f"Filter {filter_name} not found.")

    # @classmethod
    def from_modelclass(self, modelclass, **kwargs):
        SpectralElement.__init__(self, modelclass, **kwargs)

    def _database_loader(self):
        """
        Load a filter from the database.

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """
        meta, wave, flux = _read_spec(self.path, self.library.data_type,
                                      **self.library.read_kwargs)

        # TODO: why?
        wave = wave.to(u.AA)

        SpectralElement.__init__(self, Empirical1D, points=wave,
                                 lookup_table=flux.value, meta=meta)

    def _svo_loader(self, filter_name):
        wave, trans = download_svo_filter(filter_name)
        meta = {"filter_name": filter_name, "source": "SVO"}
        self.library.read_kwargs["wave_unit"] = u.Unit(wave.unit)

        SpectralElement.__init__(self, Empirical1D, points=wave,
                                 lookup_table=trans, meta=meta)

    @classmethod
    def from_vectors(
            cls,
            waves,
            trans,
            meta: dict | None = None,
            wave_unit: u.Unit = u.AA):
        """
        Create a ``Passband`` directly from from array-like.

        Parameters
        ----------
        waves: list-like
        trans: list-like
        meta: dictionary containing the metadata
        wave_unit: u.Quantity, defaulted to angstroms

        Returns
        -------
            New ``Passband`` instance.
        """
        waves <<= wave_unit
        modelclass = SpectralElement(
            Empirical1D, points=waves, lookup_table=trans, meta=meta
        )
        return cls(modelclass=modelclass)

    @classmethod
    def from_file(cls, filename: Path | str, **kwargs):
        """
        Create a ``Passband`` from a file.

        Parameters
        ----------
        filename : TYPE
            DESCRIPTION.

        Returns
        -------
            New ``Passband`` instance.
        """
        try:
            modelclass = SpectralElement.from_file(str(filename), **kwargs)
        except FileNotFoundError as err:
            raise SpextraError(f"File {filename} not found.") from err

        return cls(modelclass=modelclass)

    @classmethod
    @u.quantity_input(equivalencies=u.spectral())
    def gaussian(
            cls,
            center: u.Quantity[u.AA],
            fwhm: u.Quantity[u.AA],
            peak,
            **kwargs):
        """
        Create a ``Passband`` with a gaussian shape with given user parameters.

        Parameters
        ----------
        center : TYPE
            DESCRIPTION.
        fwhm : TYPE
            DESCRIPTION.
        peak : TYPE
            DESCRIPTION.

        Returns
        -------
            New ``Passband`` instance.
        """
        sigma = fwhm.to(u.AA) / (2.0 * np.sqrt(2.0 * np.log(2.0)))

        modelclass = SpectralElement(
            Gaussian1D,
            amplitude=peak,
            mean=center.to(u.AA),
            stddev=sigma,
            **kwargs,
        )
        return cls(modelclass=modelclass)

    @classmethod
    @u.quantity_input(equivalencies=u.spectral())
    def square(
            cls,
            wmin: u.Quantity[u.AA],
            wmax: u.Quantity[u.AA],
            transmission):
        """
        Create a ``Passband`` with a rectangular shape.

        Parameters
        ----------
        wmin : TYPE
            DESCRIPTION.
        wmax : TYPE
            DESCRIPTION.
        transmission : TYPE
            DESCRIPTION.

        Returns
        -------
            New ``Passband`` instance.
        """
        center = (wmax + wmin) * 0.5
        width = wmax - wmin
        modelclass = SpectralElement(
            Box1D,
            amplitude=transmission,
            x_0=center.to(u.AA),
            width=width.to(u.AA),
        )
        return cls(modelclass=modelclass)


class ExtinctionCurve(ReddeningLaw, ExtCurveContainer):
    """
    Class to handle extinction curves.

    TODO: Implement proper __add__, __sub__, __mul__, etc

    """

    DEFAULT_CURVES = DEFAULT_DATA.extcurves

    def __init__(self, curve_name: str | None = None, modelclass=None):
        if curve_name is not None:
            self.from_curve_name(curve_name)
        elif modelclass is not None:
            self.from_modelclass(modelclass)
        else:
            raise ConstructorError(
                "Either curve_name or modelclass must be passed."
            )

    def from_curve_name(self, curve_name: str):
        curve_name = self.DEFAULT_CURVES.get(curve_name, curve_name)

        ExtCurveContainer.__init__(self, curve_name)
        self._database_loader()

    def from_modelclass(self, modelclass):
        ReddeningLaw.__init__(self, modelclass)

    def _database_loader(self):
        """
        Load a filter from the database.

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """
        meta, wave, flux = _read_spec(self.path, self.library.data_type,
                                      **self.library.read_kwargs)

        ReddeningLaw.__init__(self, Empirical1D, points=wave,
                              lookup_table=flux, meta=meta)

    @classmethod
    def from_vectors(
            cls,
            waves,
            flux,
            meta: dict | None = None,
            wave_unit: u.Unit = u.AA):
        """
        Create a ``Passband`` directly from from array-like.

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
        waves <<= wave_unit

        modelclass = SpectralElement(
            Empirical1D, points=waves, lookup_table=flux, meta=meta
        )
        return cls(modelclass=modelclass)

    @classmethod
    def from_file(cls, filename: Path | str, **kwargs):
        modelclass = ReddeningLaw.from_file(filename, **kwargs)
        return cls(modelclass=modelclass)


class Spextrum(SourceSpectrum, SpectrumContainer):
    """
    Class to handle spectra.

    This class download, load, stores and manipulates the spectra.

    This class can be initialized with a remote file which will be downloaded
    from the database or with a ``synphot.BaseSpectrum`` or ``SourceSpectrum``.


    Parameters
    ----------
    template_name : Name of the template to download with format
    library/template e.g. "kc96/s0"
    modelclass : SourceSpectrum or BaseSpectrum

    """

    DEFAULT_SPECTRA = DEFAULT_DATA.spectra

    def __init__(self, template_name=None, modelclass=None, **kwargs):
        if template_name is not None:
            self.from_template_name(template_name, **kwargs)
        elif modelclass is not None:
            self.from_modelclass(modelclass, **kwargs)
        else:
            raise ValueError(
                "Either template_name or modelclass must be passed."
            )

    def from_template_name(self, template_name: str, **kwargs):
        template_name = self.DEFAULT_SPECTRA.get(template_name, template_name)

        SpectrumContainer.__init__(self, template_name)
        self._database_loader(**kwargs)
        self.repr = f"({self.template})"

    def from_modelclass(self, modelclass, **kwargs):
        SourceSpectrum.__init__(self, modelclass=modelclass, **kwargs)
        # FIXME: this None is problematic, as it's used for the repr, which
        #        makes an instance of this look like None o.O
        self.repr = None

    def _database_loader(self, **kwargs):
        """
        Load a template from the database.

        Returns
        -------
        meta: metadata of the spectra (header)
        lam: wavelengths
        flux: flux
        """
        meta, wave, flux = _read_spec(
            self.path, self.library.data_type, **self.library.read_kwargs
        )

        SourceSpectrum.__init__(
            self,
            Empirical1D,
            points=wave,
            lookup_table=flux,
            meta=meta,
            **kwargs,
        )

    @classmethod
    def from_arrays(
        cls,
        waves,
        flux,
        meta: dict | None = None,
        wave_unit: u.Unit = u.AA,
        flux_unit: u.Unit = units.FLAM,
    ):
        """
        Create a ``Passband`` directly from from an array-like.

        Parameters
        ----------
        waves: array-like
        flux: array-like
        meta: dictionary containing the metadata
        wave_unit: u.Quantity, defaulted to angstroms
        flux_unit: u.Quantiy, defaulted to FLAM

        Returns
        -------
        Passband
        """
        try:
            waves <<= wave_unit
        except u.UnitConversionError:
            waves = waves.to(wave_unit, equivalencies=u.spectral())
        flux <<= flux_unit

        modelclass = SourceSpectrum(
            Empirical1D, points=waves, lookup_table=flux, meta=meta
        )

        spex = cls(modelclass=modelclass)
        spex.repr = repr(spex.model)
        return spex

    @classmethod
    def from_file(cls, filename: Path | str, **kwargs):
        """
        Call ``synphot.SourceSpectrum.from_file()`` method.

        Parameters
        ----------
        filename
        kwargs

        Returns
        -------
        spex : Spextrum
            New ``Spextrum`` instance.
        """
        modelclass = SourceSpectrum.from_file(filename, **kwargs)
        spex = cls(modelclass=modelclass)
        spex.repr = f".from_file({filename})"

        return spex

    @classmethod
    @u.quantity_input
    def flat_spectrum(
        cls,
        amplitude: u.Quantity[u.ABmag] = 0 * u.ABmag,
        waves: u.Quantity[u.AA] | None = None,
    ):
        """
        Create a flat spectrum in the preferred system scaled to a magnitude.

        Defaults to a zero magnitude spectrum.

        Parameters
        ----------
        amplitude: float, u.Quantity
            amplitude/magnitude of the reference spectrum, default=0
            default is u.ABmag
            for vega use u.mag

        waves: The waveset of the reference spectrum if not Vega
           if not provided they will be created at a resolution of R~800

        Returns
        -------
        spex : Spextrum
            New ``Spextrum`` instance.
        """
        waves = (waves or _default_waves()).to(u.AA).value

        # The if-statement below also allowed amplitude.unit to be
        # u.Unit("vegamag"). Vegamag is removed from astropy, so the
        # if-statement is removed too. This assertion is added as a sanity
        # check ensure that any bugs introduced by changing the if-statement
        # is explicitly raised here, instead of silently propagated.
        assert (
            not amplitude.unit.to_string() == "vegamag"
        ), "The use of vegamag is deprecated."
        if amplitude.unit is u.Unit("mag"):
            spex = get_vega_spectrum()
            spex *= 10 ** (-0.4 * amplitude.value)
            # system_name = amplitude.unit
        else:
            const = ConstFlux1D(amplitude=amplitude)
            spex = cls(
                modelclass=Empirical1D, points=waves, lookup_table=const(waves)
            )
            # system_name = amplitude.unit

        spex.repr = f".flat_spectrum(amplitude={amplitude!r})"
        return spex

    @classmethod
    @u.quantity_input
    def black_body_spectrum(
        cls,
        temperature: u.Quantity[u.K] = 9500 * u.K,
        amplitude: u.Quantity[u.ABmag] = 0 * u.ABmag,
        filter_curve: str | Passband | SpectralElement | None = None,
        waves: u.Quantity[u.AA] | None = None,
    ):
        """
        Produce a blackbody spectrum for a given temperature.

        And scale it to a magnitude in a filter.

        Parameters
        ----------
        temperature: the temperature in Kelvin degrees
        amplitude: ``astropy.Quantity``, float
            The value that the spectrum should have in the given filter.
            Acceptable astropy quantities are:

            - u.mag : Vega magnitudes
            - u.ABmag : AB magnitudes
            - u.STmag : HST magnitudes
            - u.Jy : Jansky per filter bandpass

            Additionally the ``FLAM`` and ``FNU`` units from ``synphot.units``
            can be used when passing the quantity for `amplitude`.

        filter_curve : str
            Name of a filter from

            - a generic filter name (see ``FILTER_DEFAULTS``)
            - a spanish-vo filter service ref. (e.g. ``"Paranal/HAWKI.Ks"``)
            - a filter in the spextra database
            - a filename with the filter file
            - a ``Passband`` or ``synphot.SpectralElement`` object

        Returns
        -------
        a scaled black-body spectrum
        """
        waves = (waves or _default_waves()).to(u.AA)
        blackbody = BlackBody1D(temperature=temperature.to(u.K).value)

        spex = cls(modelclass=Empirical1D, points=waves,
                   lookup_table=blackbody(waves))
        spex = spex.scale_to_magnitude(amplitude=amplitude,
                                       filter_curve=filter_curve)
        spex.repr = (f".black_body_spectrum(amplitude={amplitude!r}, "
                     f"temperature={temperature!r}, "
                     f"filter_curve={filter_curve!r})")
        return spex

    @classmethod
    @u.quantity_input(equivalencies=u.spectral())
    def powerlaw(
        cls,
        alpha=1,
        x_0: u.Quantity[u.AA] = 5000 * u.AA,
        amplitude: u.Quantity[u.ABmag] = 0 * u.ABmag,
        filter_curve: str | Passband | SpectralElement | None = None,
        waves: u.Quantity[u.AA] | None = None,
    ):
        """
        Return a power law spectrum.

        F(lambda) ~ lambda^alpha scaled to a magnitude (amplitude) in a
        particular band.

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
            - a spanish-vo filter service ref. (e.g. ``"Paranal/HAWKI.Ks"``)
            - a filter in the spextra database
            - a filename with the filter file
            - a ``Passband`` or ``synphot.SpectralElement`` object

        Returns
        -------
        spex : Spextrum
            New ``Spextrum`` instance.
        """
        waves = (waves or _default_waves()).to(u.AA)
        power = SourceSpectrum(
            PowerLawFlux1D, amplitude=1, x_0=x_0, alpha=alpha)
        spex = cls(modelclass=Empirical1D, points=waves,
                   lookup_table=power(waves))

        spex = spex.scale_to_magnitude(
            amplitude=amplitude, filter_curve=filter_curve)

        spex.repr = (f"Spextrum.powerlaw(alpha={alpha!r}, x_0={x_0!r}, "
                     f"amplitude={amplitude!r}, filter_curve={filter_curve!r})")
        return spex

    @classmethod
    def emission_line_spectra(cls, *args, **kwargs):
        warnings.warn(
            "The 'emission_line_spectra' constructor is deprecated and will be"
            " removed in a future version. Please use the identical "
            "'emission_line_spectrum' instead.", DeprecationWarning,
            stacklevel=2)
        return cls.emission_line_spectrum(*args, **kwargs)

    @classmethod
    @u.quantity_input(equivalencies=u.spectral())
    def emission_line_spectrum(
        cls,
        center: u.Quantity[u.AA],
        fwhm: u.Quantity[u.AA],
        flux,
        amplitude: u.Quantity[u.ABmag] = 40 * u.ABmag,
        waves: u.Quantity[u.AA] | None = None,
    ):
        """
        Create a emission line spextrum superimpossed to a faint continuum.

        Parameters
        ----------
        center
        fwhm
        flux
        amplitude
        waves

        Returns
        -------
        spex : Spextrum
            New ``Spextrum`` instance.
        """
        if waves is None:
            wmin = np.min(center.to(u.AA).value) - 2000
            wmax = np.max(center.to(u.AA).value) + 2000
            step = np.min(fwhm.to(u.AA).value) / 3  # for Nyquist sampling
            waves = np.arange(wmin, wmax, step) * u.AA

        spex = cls.flat_spectrum(amplitude=amplitude, waves=waves)
        spex = spex.add_emi_lines(center=center, fwhm=fwhm, flux=flux)

        return spex

    @property
    def wave_min(self):
        return np.min(self.waveset)

    @property
    def wave_max(self):
        return np.max(self.waveset)

    def cut(self, wave_min, wave_max):
        """
        Cut the spectrum between wmin and wmax.

        Parameters
        ----------
        wave_min: float, u.Quantity,
        wave_max: float, u.Quantity,

        Returns
        -------
        spex : Spextrum
            New ``Spextrum`` instance.
        """
        wave_min <<= u.AA
        wave_max <<= u.AA

        new_waves = self.waveset[
            (self.waveset >= wave_min) & (self.waveset <= wave_max)
        ]
        spex = Spextrum(
            modelclass=Empirical1D,
            points=new_waves,
            lookup_table=self(new_waves),
            meta=self.meta,
        )
        spex = self._restore_attr(spex)

        return spex

    def redshift(self, z=0, vel=0):
        """
        Red- or blueshift a spectrum.

        Parameters
        ----------
        z : TYPE, optional
            Redshift. The default is 0.
        vel : TYPE, optional
            Radial velocity. If no unit are present it is assumed to be in m/s.
            The default is 0.

        Raises
        ------
        ValueError
            Raised if `z` or `vel` are unphysical.

        Returns
        -------
        spex : Spextrum
            New ``Spextrum`` instance.
        """
        if vel != 0:
            vel <<= u.m / u.s  # assumed to be in m/s

            z = (vel / speed_of_light).value
        if z <= -1:
            raise ValueError("Redshift or velocity unphysical")

        lam = self.waveset * (1 + z)
        flux = self(self.waveset)
        self.meta.update({"redshift": z})

        spex = Spextrum(
            modelclass=Empirical1D,
            points=lam,
            lookup_table=flux,
            meta=self.meta,
        )
        spex = self._restore_attr(spex)
        return spex

    @u.quantity_input
    def scale_to_magnitude(
        self,
        amplitude: u.Quantity[u.ABmag],
        filter_curve: str | Passband | SpectralElement | None = None,
    ):
        """
        Scale a Spectrum to a value in a filter.

        Copied from scopesim.effects.ter_curves.scale_spectrum with slight
        modifications.

        Parameters
        ----------
        amplitude : ``astropy.Quantity``, float
            The value that the spectrum should have in the given filter.
            Acceptable astropy quantities are:

            - u.mag : Vega magnitudes
            - u.ABmag : AB magnitudes
            - u.STmag : HST magnitudes
            - u.Jy : Jansky per filter bandpass

            Additionally the ``FLAM`` and ``FNU`` units from ``synphot.units``
            can be used when passing the quantity for `amplitude`:

        filter_curve : str
            Name of a filter from

            - a generic filter name (see ``DEFAULT_FILTERS``)
            - a spanish-vo filter service ref. (e.g. ``"Paranal/HAWKI.Ks"``)
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
        if Path(filter_curve).exists():
            filter_curve = Passband.from_file(filename=filter_curve)
        elif not isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = Passband(filter_curve)

        ref_spec = self.flat_spectrum(amplitude=amplitude)

        ref_flux = Observation(ref_spec,
                               filter_curve).effstim(flux_unit=units.PHOTLAM)
        real_flux = Observation(self,
                                filter_curve).effstim(flux_unit=units.PHOTLAM)

        scale_factor = ref_flux / real_flux
        spex = self * scale_factor
        spex = self._restore_attr(spex)
        spex.meta.update({"magnitude": amplitude,
                          "filter_curve": filter_curve})

        return spex

    def get_magnitude(
        self,
        filter_curve: str | Passband | SpectralElement | None = None,
        system_name: str = "AB",
    ):
        """
        Obtain the magnitude in filter for a user specified photometric system.

        Parameters
        ----------
        filter_curve : str
            Name of a filter from

            - a generic filter name (see ``DEFAULT_FILTERS``)
            - a spanish-vo filter service ref. (e.g. ``"Paranal/HAWKI.Ks"``)
            - a filter in the spextra database
            - the path to the file containing the filter (see ``Passband``)
            - a ``Passband`` or ``synphot.SpectralElement`` object

        system_name: str
            The photometric system Vega, AB or ST


        Returns
        -------
        spectrum : a Spextrum
            Input spectrum scaled to the given amplitude in the given filter.
        """
        if Path(filter_curve).exists():
            filter_curve = Passband.from_file(filename=filter_curve)
        elif not isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = Passband(filter_curve)

        if system_name.lower() in ["vega"]:
            unit = u.mag
        elif system_name.lower() in ["st", "hst"]:
            unit = u.STmag
        else:
            unit = u.ABmag

        ref_spec = self.flat_spectrum(amplitude=0 * unit)
        ref_flux = Observation(ref_spec, filter_curve).effstim(
            flux_unit=units.PHOTLAM
        )
        real_flux = Observation(self, filter_curve).effstim(
            flux_unit=units.PHOTLAM
        )

        mag = -2.5 * np.log10(real_flux.value / ref_flux.value)

        return mag * unit

    def photons_in_range(
        self,
        wmin=None,
        wmax=None,
        area=1 * u.cm**2,
        filter_curve: str | Passband | SpectralElement | None = None,
    ):
        """
        Return the number of photons between wave_min and wave_max.

        Or within a bandpass (filter).

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
            - a spanish-vo filter service ref. (e.g. ``"Paranal/HAWKI.Ks"``)
            - a filter in the spextra database
            - the path to the file containing the filter (see ``Passband``)
            - a ``Passband`` or ``synphot.SpectralElement`` object

        Returns
        -------
        counts : u.Quantity array

        """
        area <<= u.cm**2
        wmin = _angstrom_value(wmin)
        wmax = _angstrom_value(wmax)

        if filter_curve is None:
            # this makes a bandpass out of wmin and wmax
            try:
                mid_point = 0.5 * (wmin + wmax)
                width = abs(wmax - wmin)
                filter_curve = SpectralElement(
                    Box1D, amplitude=1, x_0=mid_point, width=width
                )
            except ValueError:
                raise ArgumentError("Please specify wmin/wmax or a filter")

        elif Path(filter_curve).exists():
            filter_curve = Passband.from_file(filename=filter_curve)
        elif not isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = Passband(filter_curve)

        obs = Observation(self, filter_curve)
        counts = obs.countrate(area=area)

        return counts

    def get_flux(
        self,
        wmin=None,
        wmax=None,
        filter_curve: str | Passband | SpectralElement | None = None,
        flux_unit: u.Unit = units.FLAM,
    ):
        """
        Return the flux within a passband.

        Parameters
        ----------
        wmin : float, u.Quantity
            minimum wavelength
        wmax : float, u.Quantity
            maximum wavelength
        filter_curve : str
            Name of a filter from

            - a generic filter name (see ``DEFAULT_FILTERS``)
            - a spanish-vo filter service ref. (e.g. ``"Paranal/HAWKI.Ks"``)
            - a filter in the spextra database
            - the path to the file containing the filter (see ``Passband``)
            - a ``Passband`` or ``synphot.SpectralElement`` object

        flux_unit: synphot.units, u.Quantity

        Returns
        -------
        flux

        """
        wmin = _angstrom_value(wmin)
        wmax = _angstrom_value(wmax)

        if filter_curve is None:
            # this makes a bandpass out of wmin and wmax
            try:
                filter_curve = Passband.square(wmin=wmin, wmax=wmax)
            except ValueError:
                raise ArgumentError("Please specify wmin/wmax or a filter")
        elif Path(filter_curve).exists():
            filter_curve = Passband.from_file(filename=filter_curve)
        elif not isinstance(filter_curve, (Passband, SpectralElement)):
            filter_curve = Passband(filter_curve)

        flux = Observation(self, filter_curve).effstim(flux_unit=flux_unit)

        return flux

    def add_emi_lines(self, center, fwhm, flux):
        """
        Add emission lines to an `Spextrum`.

        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)

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
        # Converting flux is done by GaussianFlux1D anyway...
        centers = list(always_iterable(_angstrom_value(center)))
        fluxes = list(always_iterable(flux))
        fwhms = list(always_iterable(_angstrom_value(fwhm)))

        self.meta.update({"em_lines": {"center": centers,
                                       "flux": fluxes,
                                       "fwhm": fwhms}})

        for cen, flx, fwh in zip(centers, fluxes, fwhms):
            line = GaussianFlux1D(mean=cen, total_flux=flx, fwhm=fwh)
            # lam = line.sampleset(factor_step=0.3)  # bit better than Nyquist
            lam = self.waveset.to(u.AA).value
            g_em = SourceSpectrum(Empirical1D, points=lam,
                                  lookup_table=line(lam))
            self += g_em

        spex = self._restore_attr(Spextrum(modelclass=self))
        return spex

    def add_abs_lines(self, center, ew, fwhm):
        """
        Add an absorption line to a spectrum.

        With center, fwhm and equivalent width specified by the user.
        It also supports emission lines if `ew` is negative.

        TODO: accept different profiles (Lorentz1D, Voigt1D, etc)

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
        centers = list(always_iterable(_angstrom_value(center)))
        eqws = list(always_iterable(_angstrom_value(ew)))
        fwhms = list(always_iterable(_angstrom_value(fwhm)))

        self.meta.update({"em_lines": {"center": centers,
                                       "ew": eqws,
                                       "fwhm": fwhms}})

        for cent, eqw, fwh in zip(centers, eqws, fwhms):
            # to keep the convention that EL are negative and ABS are positive
            sign = -1 * np.sign(eqw)
            left, right = cent - np.abs(eqw / 2), cent + np.abs(eqw / 2)
            wavelengths = self.waveset[(self.waveset.value >= left) &
                                       (self.waveset.value <= right)]
            fluxes = units.convert_flux(wavelengths=wavelengths,
                                        fluxes=self(wavelengths),
                                        out_flux_unit=units.FLAM)
            flux = np.trapz(fluxes.value, wavelengths.value)
            line = GaussianFlux1D(mean=cent, total_flux=sign * flux,
                                  fwhm=fwh)
            lam = line.sampleset(factor_step=0.35)  # bit better than Nyquist
            flux = line(lam).clip(min=0.0)
            g_abs = SourceSpectrum(Empirical1D, points=lam,
                                   lookup_table=flux)

            self += g_abs

            if (self(wavelengths).value < 0).any():
                warnings.warn("Warning: Flux<0 for specified EW and FHWM, "
                              "setting it to Zero")
                waves = self.waveset[self(self.waveset) < 0]
                zero_sp = SourceSpectrum(Empirical1D, points=waves,
                                         lookup_table=-1 * self(waves).value)
                self += zero_sp

        spex = self._restore_attr(Spextrum(modelclass=self))
        return spex

    @staticmethod
    def _get_extinct(curve_name: str, Ebv, Av, Rv):
        if Av is not None:
            Ebv = Av / Rv
        extcurve = ExtinctionCurve(curve_name)
        return extcurve.extinction_curve(Ebv)

    def redden(
        self,
        curve_name: str,
        Ebv: float = 0.0,
        Av: float | None = None,
        Rv: float = 3.1,
    ):
        """
        Attenuate a spectrum with a extinction curve normalized to a E(B-V).

        Parameters
        ----------
        curve_name: str, name of the extinction curve
        Av: float, Av parameter
        Ebv: float, E(B-V) color excess

        Returns
        -------
        spex : Spextrum
            New attenuated ``Spextrum`` instance.
        """
        extinction = self._get_extinct(curve_name, Ebv, Av, Rv)
        spex = Spextrum(modelclass=self.model * extinction.model)
        spex = self._restore_attr(spex)
        spex.meta.update({"Reddening with": curve_name,
                          "Reddening with E(B-V)": Ebv})
        return spex

    def deredden(
        self,
        curve_name: str,
        Av: float | None = None,
        Ebv: float = 0.0,
        Rv: float = 3.1,
    ):
        """
        De-redden a spectrum.

        Parameters
        ----------
        curve_name: str, name of the extinction curve
        Av: float, Av parameter
        Ebv: float, E(B-V) color excess

        Returns
        -------
        spex : Spextrum
            New attenuated ``Spextrum`` instance.
        """
        extinction = self._get_extinct(curve_name, Ebv, Av, Rv)
        spex = Spextrum(modelclass=self.model / extinction.model)
        spex = self._restore_attr(spex)
        spex.meta.update({"Dereddening with": curve_name,
                          "Dereddening with E(B-V)": Ebv})
        return spex

    @u.quantity_input(equivalencies=u.spectral())
    def rebin_spectra(self, new_waves: u.Quantity[u.AA]):
        """
        Rebin a synphot spectra to a new wavelength grid conserving flux.

        Grid does not need to be linear and can be at higher or lower
        resolution.

        Parameters
        ----------
        new_waves: an array of the output wavelenghts in Angstroms but other
        units can be specified.

        Returns
        -------
        spex : Spextrum
            New rebinned ``Spextrum`` instance.
        """
        new_waves = new_waves.to(u.AA).value

        waves = self.waveset.value  # else assumed to be angstroms
        flux = np.ones_like(waves)
        filt = SpectralElement(Empirical1D, points=waves, lookup_table=flux)
        obs = Observation(self, filt, binset=new_waves, force="taper")
        newflux = obs.binflux
        rebin_spec = Empirical1D(
            points=new_waves, lookup_table=newflux, meta=self.meta
        )
        spex = Spextrum(modelclass=rebin_spec)
        spex = self._restore_attr(spex)
        spex.meta.update({"rebinned_spectra": True})

        return spex

    def logrebin(self):
        """
        TODO: Estimate optimal rebinning factors here.
        Returns
        -------
        a Spextrum with a log rebinned in wavelength
        """
        nbins = int(self.waveset.ptp() / np.diff(self.waveset).min())
        logwaves = np.geomspace(self.wave_min, self.wave_max, nbins)
        return self.rebin_spectra(logwaves)

    def smooth(self, sigma):
        """
        Smooth the Spectrum with a Gaussian Kernel.

        Expressed in km/s (or other velocity unit).

        The spectra is first log-rebinned so it has a constant velocity spacing
        and after that is smoothed with a gaussian kernel.


        Parameters
        ----------
        sigma: u.Quantity, width of the Kernel in

        Returns
        -------
        smoothed spextrum
        """
        sigma <<= u.km / u.s

        sp_log = self.logrebin()
        step_size = np.median(np.diff(sp_log.waveset) / sp_log.waveset[:1]
                              * speed_of_light).to(u.km/u.s)

        if step_size > sigma:
            warnings.warn("Spectrum is undersampled for provided sigma value, "
                          f"{step_size} > {sigma}.")

        conv_sigma = sigma / step_size
        smoothed_flux = gaussian_filter1d(
            sp_log(sp_log.waveset), conv_sigma.value)

        self.meta.update({"KERNEL_SIZE": sigma.value})
        smoothed = Empirical1D(
            points=sp_log.waveset.value, lookup_table=smoothed_flux,
            meta=self.meta)
        spex = self._restore_attr(Spextrum(modelclass=smoothed))
        return spex

    @u.quantity_input(equivalencies=u.spectral())
    def add_noise(self, wmin: u.Quantity[u.AA], wmax: u.Quantity[u.AA]):
        """Return spectrum with given S/N in the specified wavelength range."""
        raise NotImplementedError()

    def _restore_attr(self, spextrum):
        """
        just a hack to down stream important attributes in particular methods.
        there must be a better way
        """
        temp_dict = {
            key: value
            for key, value in self.__dict__.items()
            if not key.startswith("_")
        }
        spextrum.__dict__.update(temp_dict)
        return spextrum

    # ------ Copied from synphot.SourceSpectrum so operations can also happen here --------
    # docstrings modified...

    def __add__(self, other):
        """Add `self` with `other`."""
        self._validate_other_add_sub(other)
        result = self.__class__(modelclass=self.model + other.model)
        self._merge_meta(self, other, result)
        return result

    def __sub__(self, other):
        """Subtract `other` from `self`."""
        self._validate_other_add_sub(other)
        result = self.__class__(modelclass=self.model - other.model)
        self._merge_meta(self, other, result)
        return result

    def __mul__(self, other):
        """Multiply `self` and `other`."""
        self._validate_other_mul_div(other)

        if isinstance(other, (u.Quantity, numbers.Number)):
            newcls = self.__class__(modelclass=self.model | Scale(other))
        elif isinstance(other, BaseUnitlessSpectrum):
            newcls = self.__class__(modelclass=self.model * other.model)
        else:  # Source spectrum
            raise exceptions.IncompatibleSources(
                "Cannot multiply two source spectra together."
            )

        self._merge_meta(self, other, newcls)
        return newcls

    def _validate_other_add_sub(self, other):
        """
        Check conditions for other to satisfy before add/sub.

        Copied from ``synphot.SourceSpectrum`` so additions can happen with a
        ``SourceSpectrum``.
        """
        if not isinstance(other, (self.__class__, SourceSpectrum)):
            raise exceptions.IncompatibleSources(
                f"Can only operate on {self.__class__.__name__}."
            )

    def _validate_other_mul_div(self, other):
        """
        Check conditions for other to satisfy before mul/div.

        Copied from ``synphot.SourceSpectrum`` so multiplications can happen
        with a ``SourceSpectrum``.
        """

        if not isinstance(
            other,
            (
                u.Quantity,
                numbers.Number,
                BaseUnitlessSpectrum,
                SourceSpectrum,
                self.__class__,
            ),
        ):
            raise exceptions.IncompatibleSources(
                "Can only operate on scalar number/Quantity or spectrum."
            )
        if isinstance(other, u.Quantity) and (
            other.unit.decompose() != u.dimensionless_unscaled
            or not np.isscalar(other.value)
            or not isinstance(other.value, numbers.Real)
        ):
            raise exceptions.IncompatibleSources(
                "Can only operate on real scalar dimensionless Quantity."
            )
        if isinstance(other, numbers.Number) and not (
            np.isscalar(other) and isinstance(other, numbers.Real)
        ):
            raise exceptions.IncompatibleSources(
                "Can only operate on real scalar number."
            )

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"{self.__class__.__name__}{self.repr}"


# ------------------------------ END    -------------------------------------------


def get_vega_spectrum():
    """
    Retrieve the Vega spectrum from the database.

    Notes
    -----
    To access wavelength and fluxes use::
        spex = get_vega_spectrum()
        waves, fluxes = spex.waveset, spex(spex.waveset)

    """
    return Spextrum("ref/vega")


def _default_waves():
    """Generate default waveset with R~805."""
    waves, _ = utils.generate_wavelengths(
        minwave=100, maxwave=50000, num=5000, log=True, wave_unit=u.AA
    )
    return waves


def _read_spec(fname, fmt, **kwargs):
    if fmt == "fits":
        return read_fits_spec(fname, ext=1, **kwargs)
    elif fmt == "ascii":
        kwargs = {
            key: kwargs[key]
            for key in kwargs.keys() & {"wave_unit", "flux_unit"}
        }
        return read_ascii_spec(fname, **kwargs)
    else:
        raise SpextraError("invalid data type")


def _angstrom_value(value):
    if isinstance(value, u.Quantity):
        return value.to(u.AA, equivalencies=u.spectral()).value
    return value
