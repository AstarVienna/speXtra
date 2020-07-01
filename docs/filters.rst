.. _filters

*******
Filters
*******


At the moment ``speXtra`` mostly relies on `tynt <https://tynt.readthedocs.io/en/latest/>` to fetch
the filters from the `SVO database <http://svo2.cab.inta-csic.es/theory/fps/>`, however a few
filters systems have been included in the database

To list the filter systems included in the database

.. code-block:: python

    from spextra import SpecDatabase
    db = SpecDatabase()
    print(db.filter_systems)

To list the filter systems available in the `SVO database <http://svo2.cab.inta-csic.es/theory/fps/>`

.. code-block:: python

    from spextra.database import get_filter_systems
    print(get_filter_systems())

--------------------------------------------------------------------

.. _micado

MICADO
======

.. literalinclude:: ../database/filter_systems/micado/index.yml
    :language: yaml

