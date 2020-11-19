.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
      th { display:none; }
    </style>


*******
speXtra
*******


.. image:: https://travis-ci.org/miguelverdugo/speXtra.svg?branch=master
    :target: https://travis-ci.org/github/miguelverdugo/speXtra
    :alt: Tests Status

.. image:: https://readthedocs.org/projects/spextra/badge/?version=latest
    :target: https://spextra.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

speXtra
=======

``speXtra`` is a :python:`python <>` library for managing and manipulating astronomical spectra.

Under the hood it packages many :synphot:`synphot <>` workflows to perform typical
treatment of astronomical spectra, including scaling, binning, smoothing, redshifting, etc.
``speXtra`` does not make measurement on the spectra, except obtaining magnitudes.

``speXtra`` comes with a built-in :ref:`database <database>` to facilitate the access
to commonly used spectral templates, extinction curves and filters.


.. toctree::
   :maxdepth: 1


   install
   starting
   database
   spectral_libraries
   extinction_curves
   filters










:doc:`reference`
================

More...
=======

.. toctree::
   :maxdepth: 1
   :titlesonly:


   about

