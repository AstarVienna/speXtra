# -*- coding: utf-8 -*-
import os
import inspect
from . import spextra
from . import database

from .spextra import Spextrum, make_passband
from .database import SpecDatabase

__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))
__data_dir__ = os.path.join(__pkg_dir__, "data")

