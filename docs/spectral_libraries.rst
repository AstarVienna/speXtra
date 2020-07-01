.. _spectral-libraries

******************
Spectral Libraries
******************



A number of spectral libraries have been include in ``speXtra``.
See :ref:`database-contents` for a summary what is included in the database.

To list the names of the libraries included in the database

.. code-block:: python

    from spextra import SpecDatabase
    db = SpecDatabase()
    print(db.libraries)


To see which templates are available in each library

.. code-block:: python

    from spextra.database import SpecLibrary
    lib = SpecLibrary(name)
    print(lib.templates)


Below you can find a detailed description of each library.

--------------------------------------------------------------------

.. _kc96

The Kinney-Calzetti Spectral Atlas of Galaxies
==============================================

.. literalinclude:: ../database/libraries/kc96/index.yml
    :language: yaml

.. _pickles

Pickles Stellar Library
=======================

.. literalinclude:: ../database/libraries/pickles/index.yml
    :language: yaml


SDSS galaxy composite spectra
=============================

.. literalinclude:: ../database/libraries/sdss/index.yml
    :language: yaml


IRTF spectral library
=====================

.. literalinclude:: ../database/libraries/irtf/index.yml
    :language: yaml


AGN templates
=============

.. literalinclude:: ../database/libraries/agn/index.yml
    :language: yaml


Emission line nebulae
=====================


.. literalinclude:: ../database/libraries/nebulae/index.yml
    :language: yaml

Galaxy SEDs from the UV to the Mid-IR
=====================================

.. literalinclude:: ../database/libraries/brown/index.yml
    :language: yaml


Subset of Kurucz 1993 Models
============================

.. literalinclude:: ../database/libraries/kurucz/index.yml
    :language: yaml


Supernova Legacy Survey
=======================

.. literalinclude:: ../database/libraries/sne/index.yml
    :language: yaml


Flux/Telluric standards with X-Shooter
=======================================

.. literalinclude:: ../database/libraries/moehler/index.yml
    :language: yaml