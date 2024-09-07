# -*- coding: utf-8 -*-
from . import spextra
from . import database

from .spextra import Spextrum, Passband, ExtinctionCurve, get_vega_spectrum
from .database import spextra_database, DEFAULT_DATA
from .libraries import SpecLibrary, FilterSystem, ExtCurvesLibrary
from .configuration import config
from .exceptions import SpextraError
