.. _database:


********
Database
********


``speXtra`` comes with a built-in database of spectral templates, extinction curves and astronomical filters.

The database organized through ``yaml`` files which describe the contents of the different data files.

The inner workings of the database are transparent to the user and generally the user does not
need to deal with the database when working with ``speXtra``.

There are however few things that might be useful to the user, specially when working in interactive mode
(e.g. Jupyter Notebooks)

The database is organized within a directory tree. Depending on what you are requesting, a spectra template,
a extinction curve or a filter, it will look for it at that particular place. The syntax is standard:

- ``"library_name/template_name"`` for a spectral template
- ``"extinction_curve_family/extinction_curve_name"`` for extinction curves
- ``"filter_system/filter_name"`` for astronomical filters

Below you can find the contents of the database :ref:`database-contents`


Browsing the database
=====================

There are few ways to interact with the database to examine its contents or to use them
programately








.. code-block:: python

    from spextra.database import SpecDatabase
    db = SpecDatabase()
    print(db.libraries)

it will list the spectral libraries available

.. code-block:: python

    print(db.extinction_curves)
    print(db.filter_systems)

will print the extinction curves and filter systems available


:class:`SpecLibrary` is also important to examine the contents
of each spectral library. To use it, simply call it like this

.. code-block:: python
    from spextra.database import SpecLibrary
    lib = SpecLibrary("kc96)
    print(lib.templates)

and it will list all templates available for that library.


Similarly :class:`FilterSystem` holds the information for a particular filter system and :class:`ExtinctionCurveLibrary`
for a particular extinction curve family.


--------------------------------------------------------------------

.. _database-contents:

Database contents
=================


.. literalinclude:: ../database/index.yml
    :language: yaml







.. toctree::
   :maxdepth: 2


   spectral_libraries
   extinction_curves
   filters




