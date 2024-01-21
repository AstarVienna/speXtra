.. _install:

************
Installation
************

The preferred method to install  ``spextra`` is using ``pip``

.. code-block:: bash

    pip install spextra


To install the development version, just clone the repository


.. code-block:: bash

    git clone https://github.com/AstarVienna/speXtra.git
    cd speXtra
    pip install -e .




Dependencies
------------

The following dependencies are necessary to run ``speXtra``

- `numpy <http://www.numpy.org/>`_
- `astropy <http://www.astropy.org>`_
- `scipy <http://www.scipy.org/>`_
- `synphot <http://synphot.readthedocs.io>`_
- `PyYAML <https://pyyaml.org/>`_

Additionally, you may need the following libraries for specific purposes.

- `matplotlib <http://www.matplotlib.org/>`_ for plotting
- `specutils <specutils.readthedocs.io/>`_ optional for loading external spectra in other formats
