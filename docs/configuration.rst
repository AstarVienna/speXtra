**************
Data Directory
**************


``speXtra`` will create and use an ``spextra`` subdirectory in
the AstroPy cache directory for this purpose. For example,
``$HOME/.astropy/cache/spextra``. This directory structure will match
the database in the server. However ``speXtra`` will download
the data on demand.


Changing the Data Directory
===========================

Sometimes might be necessary to change the data directory for better accessibility.
This can be done by using the following commands

.. code-block:: python

   from spextra.utils import Conf

   conf = Conf(datadir="path/to/new/datadir")
   print(conf.get_data_dir())

   # To set it back to the default

   Conf(datadir=conf.default_data_dir)


