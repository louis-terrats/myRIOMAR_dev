"""
Detects river plumes and generates associated plots.
"""

from .main_functions import (
    apply_plume_mask,
    make_and_plot_time_series_of_plume_areas
)

__all__ = [
    "apply_plume_mask",
    "make_and_plot_time_series_of_plume_areas"
]

