# myRIOMAR/CODES/__init__.py

"""
Main package initializer for myRIOMAR.
Exposes the core submodules.
"""

from . import _0_data_downloading
from . import _1_data_validation
from . import _2_regional_maps
from . import _3_plume_detection
from . import _4_X11_analysis
from . import _5_Figures_for_article
from . import TASKS
from . import DATA
      
__all__ = [
    "_0_data_downloading",
    "_1_data_validation",
    "_2_regional_maps",
    "_3_plume_detection",
    "_4_X11_analysis",
    "_5_Figures_for_article",
    "TASKS",
    "DATA",
]
