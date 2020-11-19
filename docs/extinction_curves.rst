.. _extinction_curves:

*****************
Extinction Curves
*****************



Few extinction curves have been included in ``speXtra``


To list the names of the extinction curves included in the database

.. jupyter-execute::

    from spextra import Database
    db = Database()
    db.extinction_curves

--------------------------------------------------------------------

.. _gordon:

Gordon
======

.. literalinclude:: ../database/extinction_curves/gordon/index.yml
    :language: yaml

.. _cardelli:

Cardelli
========

.. literalinclude:: ../database/extinction_curves/cardelli/index.yml
    :language: yaml


.. _calzetti:

Calzetti
========

.. literalinclude:: ../database/extinction_curves/calzetti/index.yml
    :language: yaml

