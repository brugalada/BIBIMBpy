"""
Backwards Integration Basic Interface Module for Bars with Python

This is a simple interface with AGAMA intended to simplify the process of setting up and running backward integrations. The 

BIBIMBpy requires AGAMA, numpy and scipy and matplotlib.
"""

__version__ = "1.0"

try:
    import agama
except ImportError:
    raise ImportError('AGAMA does not seem to be installed.')

try:
    import numpy
except ImportError:
    raise ImportError('NumPy does not seem to be installed.')

try:
    import scipy
except ImportError:
    raise ImportError('SciPy does not seem to be installed.')

try:
    import matplotlib
except ImportError:
    raise ImportError('MatPlotLib does not seem to be installed.')

from . import orbits
from . import initialize
from . import utils
