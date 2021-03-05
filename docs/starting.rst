.. _start:

***************
Getting Started
***************

The core of ``speXtra`` is the :class:`Spextrum` which is a wrapper of the :class:`synphot.SourceSpectrum`
with added functionalities.

For example to load a `S0` galaxy templates from the  :ref:`kc96` (`kc96`) just type

.. jupyter-execute::

    from spextra import Spextrum
    sp = Spextrum("kc96/s0")
    sp.plot()


The last statement will create a plot of the spectrum for a quick visualization



All operations available in   :class:`synphot.SourceSpectrum` are possible:

.. jupyter-execute::

    sp1 = Spextrum("kc96/s0")
    sp2 = Spextrum("agn/qso")
    sp = sp1 + 0.3*sp2
    sp.plot()



Adding emission lines
----------------------

It is possible to add emission lines, either individually or as a list. Parameters are center, flux and fwhm
`astropy.Units` are allowed. If units are not specified it defaults to Angstroms for wavelengths (center and fwhm)
and ergs/s/AA/cm^2 (FLAM) for flux.

.. jupyter-execute::

    import astropy.units as u

    sp3 = sp1.add_emi_lines(center=4000,flux=4e-13, fwhm=5*u.AA)
    sp3.plot()



Scaling to a magnitude
----------------------

.. jupyter-execute::

    sp1 = Spextrum("kc96/s0")
    sp2 = sp1.scale_to_magnitude(amplitude=13 * u.ABmag, filter_curve="g")
    sp2.plot()



Obtaining magnitudes from spectra
---------------------------------

.. jupyter-execute::

    sp1.get_magnitude(filter_curve="g")


Redshifting the spectra
------------------------

It is possible to specify a redshift `z` ir a velocity `vel` (negative velocities are allowed)

.. jupyter-execute::

    sp3 = sp2.redshift(z=1)

    import astropy.units as u

    vel = -1000 * u.km / u.s
    sp2 = sp1.redshift(vel=vel)




Smooth the spectral
-------------------

Spectra can be smoothed with a kernel with a size in velocities (default km/s)

.. jupyter-execute::

    sp1 = Spextrum("nebulae/pn")

    sigma = 500*(u.km / u.s)
    sp2 = sp1.smooth(sigma=sigma)













