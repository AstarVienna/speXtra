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


.. warning:: July 2022: The downloadable content server was retired and the data migrated to a new server.

   SpeXtra v0.25 and above have been redirected to a new server URL.

   For older versions, please either upgrade to the latest version (``pip install --upgrade spextra``), or follow these `instructions to update the server URL <https://astarvienna.github.io/server_upgrade_instructions.html>`_ in the config file.


.. toctree::
   :maxdepth: 1


   install
   starting
   notebooks/Database.ipynb
   spectral_libraries
   extinction_curves
   filters
   package






   









:doc:`reference`
================

More...
=======

.. toctree::
   :maxdepth: 1
   :titlesonly:


   about

