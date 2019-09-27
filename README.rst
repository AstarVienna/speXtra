=======
Spyctra
=======

A python tool to manage and manipulate astronomical spectra



Description
===========

Spyctra is a python tool to download, load, display and manipulate spectra of astronomical sources. It has developed to provide spectral sources to SpecSim but it may be helpful for other purposes too.

It heavily piggybacks in funcionalities provided by synphot and astropy

Functionalities
===============

- Download spectra from known sources and return it in synphot format
- Load spectra from files in different formats and return it in synphot format
- Manipulate spectra:
   - rebinning (interpolate), including log-rebin and constant velocity
   - smooth with a kernel, including a wavelengths changing kernel (think about it later)
   - redshifting (blueshifting)
   - attenuating
   - scaling to a magnitude?
- plot the spectra
- save it the wcs1d format




