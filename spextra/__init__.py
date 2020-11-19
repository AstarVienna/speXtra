# -*- coding: utf-8 -*-
import os
import inspect
from . import spextra
from . import database
from . import utils

from .spextra import Spextrum, Passband, ExtinctionCurve, DEFAULT_FILTERS, DEFAULT_SPECTRA, DEFAULT_CURVES
from .database import Database, SpecLibrary, FilterSystem, ExtCurvesLibrary, DefaultData

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

