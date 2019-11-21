# -*- coding: utf-8 -*-
"""
    Setup file for spextra.
"""
import sys
import setuptools
import pytest


def setup_sp():
    setuptools.setup(
        name="speXtra",
        description="Tool to  manage and manipulate  astronomical spectra",
        author="Miguel Verdugo",
        author_email="mverduol@gmail.com",
        license="MIT",
        url="https://github.com/miguelverdugo/speXtra",
        package_dir={'spextra': 'spextra'},
        packages=['spextra'],
        package_data={'spextra': ['spextra/data/*']},
        include_package_data=True,
        install_requires=["numpy",
                          "astropy",
                          "synphot>=0.2",
                          "matplotlib>1.5.0",
                          "tynt",
                          "pyyaml", ],
        classifiers=["Programming Language :: Python :: 3",
                     "License :: OSI Approved :: MIT License",
                     "Operating System :: OS Independent",
                     "Intended Audience :: Science/Research",
                     "Topic :: Scientific/Engineering :: Astronomy", ]
    )


if __name__ == "__main__":
    setup_sp()
