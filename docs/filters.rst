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


Many filters in the SVO database have been defined using shortcuts for quicker access.
Here is a list.

.. jupyter-execute::

    from spextra import DEFAULT_FILTERS
    DEFAULT_FILTERS


These filters can be used simply using the shortcut

.. jupyter-execute::

    from spextra import Passband

    filt = Passband("g")  # for sdss g-band filter
    filt.plot()


Additionally below you can find the filter systems for prospective ELT instruments

--------------------------------------------------------------------

.. _micado:

MICADO
======

.. literalinclude:: ../database/filter_systems/elt/micado/index.yml
    :language: yaml

METIS
=====

.. literalinclude:: ../database/filter_systems/elt/metis/index.yml
    :language: yaml
