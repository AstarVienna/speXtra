# -*- coding: utf-8 -*-
import os
import inspect
from . import spextra
from . import database
from . import utils

from .spextra import Spextrum, Passband, ExtinctionCurve
from .database import Database, SpecLibrary, FilterSystem, ExtCurvesLibrary

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

Database(reload=True, silent=True)  # ensures that the database is always up to the date


