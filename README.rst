=======
SpeXtra
=======

A python tool to manage and manipulate astronomical spectra



Description
===========

speXtra is a python tool to download, load, display and manipulate spectra of astronomical sources.
It has developed to provide spectral sources to SpecSim but it may be helpful for other purposes too.

speXtra stands in the shoulder of giants: ``synphot`` and ``astropy``


Functionalities
===============

speXtra is able to:

- Download spectra from known sources and return it in format compatible with ``synphot`` format


- Load spectra from files in different formats and return it in ``synphot`` format
   - Ascii and fits tables are natively supported by ``synphot``
   - fits files and other type of spectra will be supported by converting
     a ``specutils.Spectrum1D`` object to a ``synphot.SourceSpectrum``. Leave it to
     ``specutils`` to try to do the file loading.

- Manipulate spectra:
   - rebinning (interpolate), including log-rebin and constant velocity
   - smooth with a kernel, including a wavelength changing kernel (think about it later)
   - redshifting (blueshifting)
   - (de-)reddening and (de-)attenuating
   - scaling to a magnitude in filter
   - obtaining magnitude in filter
   - Add emission/absorption line (currently only gaussian profile is supported)

- plot the spectra



TODO
====

- Create a logo ;-)


