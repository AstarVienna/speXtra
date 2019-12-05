=======
speXtra
=======

This is the documentation of **speXtra**.

speXtra is a python module to manage and manipulate astronomical spectra. It stands on
the shoulders of giants, `astropy`_ and `synphot`_. In particular, it includes
many `synphot`_ algorithms and workflows to make easier the manipulation of spectra.

.. note::

    **speXtra** is being developed with the main purpose to become the spectral engine
    for future telescope simulation software like `ScopeSim`_. However it may be
    useful for many other purposes.

    **speXtra** is compatible only with python 3.5 and above


Installation
============

At the moment, the best way to install it is to download or clone the package, e.g.

    wget https://github.com/miguelverdugo/speXtra/archive/master.zip

    unzip master.zip

    python3 setup.py install

or

    pip install git+https://github.com/miguelverdugo/speXtra.git

In the near future, there will also be a PyPi version

The following dependencies are necessary to run **speXtra**:

    - astropy
    - synphot
    - numpy
    - PyYAML
    - tynt
    - matplotlib (if you want to plot the spectra)
    - specutils (optional)


Contents
========

.. toctree::
   :maxdepth: 2

   License <license>
   Authors <authors>
   Changelog <changelog>
   Module Reference <api/modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _toctree: http://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html
.. _reStructuredText: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _references: http://www.sphinx-doc.org/en/stable/markup/inline.html
.. _Python domain syntax: http://sphinx-doc.org/domains.html#the-python-domain
.. _astropy:  https://docs.astropy.org/en/stable/
.. _synphot: https://synphot.readthedocs.io/en/latest/
.. _ScopeSim: https://scopesim.readthedocs.io/en/latest/