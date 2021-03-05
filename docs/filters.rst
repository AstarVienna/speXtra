.. _filters:

*******
Filters
*******

Most of the filters available to ``speXtra`` are downloaded from
the `SVO database <http://svo2.cab.inta-csic.es/theory/fps/>`. However, many other
system are not included in that database so a few (at the moment) are included in the database.

To list the filter systems included in the database

.. jupyter-execute::

    from spextra import Database
    db = Database()
    db.filter_systems

The ::class::`Pasband` holds the information of an astronomical filter and can be called simply using:

.. jupyter-execute::

    from spextra import Passband
    filt = Passband("elt/micado/Y")
    filt.plot()


Many standard filters have been defined using shortcuts for quicker access.
Here is a list.


.. jupyter-execute::

    from spextra import DEFAULT_FILTERS
    DEFAULT_FILTERS


These filters can be used simply using the shortcut

.. jupyter-execute::

    from spextra import Passband

    filt = Passband("g")  # for sdss g-band filter
    filt.plot()


All methods in ::class::`Spextrum` that require to provide a filter name, can be used in either way,
using the full name or just a shortcut.

For a full list of the filters available to ``spextra`` please visit the `SVO <http://svo2.cab.inta-csic.es/theory/fps/>`_


Additionally below you can find the filter systems for prospective ELT instruments and others not
available at the SVO

--------------------------------------------------------------------

.. _micado:

MICADO
======

.. literalinclude:: ../database/filter_systems/elt/micado/index.yml
    :language: yaml


.. _metis:

METIS
=====

.. literalinclude:: ../database/filter_systems/elt/metis/index.yml
    :language: yaml

----------------------------------------------------------------------

Also the standard filters used by the ESO Exposure Time Calculator

.. _etc:

ETC
===

.. literalinclude:: ../database/filter_systems/etc/index.yml
    :language: yaml

