.. _spextra-database:


********
Database
********


.. _database:

Database
========

Create a spectra using the built-in ``'s0'`` galaxy source from the Kinney-Calzetti library
``'kc96'``:

    >>> from speXtra import Spextrum
    >>> sp = Spextrum("kc96/s0")

and ``spextra`` will download the spectra and create an object

--------------------------------------------------------------------

.. _database-contents:

Database contents
=================


.. literalinclude:: ../database/index.yml
    :language: yaml



.. _spectral-libraries

--------------------------------------------------------------------

The Kinney-Calzetti Spectral Atlas of Galaxies
==============================================


.. literalinclude:: ../database/libraries/kc96/index.yml
    :language: yaml

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