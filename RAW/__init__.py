from __future__ import print_function, division

# dirty hack solution for pickle compatibility problem
import os
import sys

RAW_DIR = os.path.dirname(os.path.abspath(__file__))
if RAW_DIR not in sys.path:
    sys.path.append(RAW_DIR)
from .RAWWrapper import RAWSimulator
from .RAWAnalysisWrapper import RAWAnalysisSimulator
from . import SASM

# TODO: Fix absolute import between python 2/3 and pickle

__NAME__ = 'RAW'
__VERSION__ = '1.2.2'
