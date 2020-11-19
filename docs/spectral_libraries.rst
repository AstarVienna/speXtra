.. _spectral-libraries:

******************
Spectral Libraries
******************



A number of spectral libraries have been include in ``speXtra``.
See :ref:`database-contents` for a summary what is included in the database.

To list the names of the libraries included in the database

.. jupyter-execute::

    from spextra import Database
    db = Database()
    db.libraries


To see which templates are available in each library

.. jupyter-execute::

    from spextra.database import SpecLibrary
    name = "kc96"
    lib = SpecLibrary(name)
    lib.templates


Below you can find a detailed description of each library.

--------------------------------------------------------------------

.. _ref:

A Library of Reference Stars
----------------------------

.. literalinclude:: ../database/libraries/ref/index.yml
    :language: yaml



.. _kc96:

The Kinney-Calzetti Spectral Atlas of Galaxies
----------------------------------------------

.. literalinclude:: ../database/libraries/kc96/index.yml
    :language: yaml

.. _pickles:

Pickles Stellar Library
-----------------------

.. literalinclude:: ../database/libraries/pickles/index.yml
    :language: yaml


.. _sdss:

SDSS galaxy composite spectra
-----------------------------

.. literalinclude:: ../database/libraries/dobos/index.yml
    :language: yaml

.. _irtf:

IRTF spectral library
---------------------

.. literalinclude:: ../database/libraries/irtf/index.yml
    :language: yaml

.. _agn:

AGN templates
-------------

.. literalinclude:: ../database/libraries/agn/index.yml
    :language: yaml

.. _nebulae:

Emission line nebulae
---------------------


.. literalinclude:: ../database/libraries/nebulae/index.yml
    :language: yaml

.. _brown:

Galaxy SEDs from the UV to the Mid-IR
-------------------------------------

.. literalinclude:: ../database/libraries/brown/index.yml
    :language: yaml


.. _kurucz:

Subset of Kurucz 1993 Models (subset)
-------------------------------------

.. literalinclude:: ../database/libraries/kurucz/index.yml
    :language: yaml


.. _sne:

Supernova Legacy Survey
-----------------------

.. literalinclude:: ../database/libraries/sne/index.yml
    :language: yaml


.. _moehler:

Flux/Telluric standards with X-Shooter
--------------------------------------

.. literalinclude:: ../database/libraries/moehler/index.yml
    :language: yaml


.. _madden:

High-Resolution Spectra of Habitable Zone Planets (example)
-----------------------------------------------------------

.. literalinclude:: ../database/libraries/madden/index.yml
    :language: yaml
