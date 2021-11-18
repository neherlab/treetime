from __future__ import print_function, division, absolute_import
version="0.8.5"

class TreeTimeError(Exception):
    """TreeTimeError class"""
    pass

class MissingDataError(TreeTimeError):
    """MissingDataError class raised when tree or alignment are missing"""
    pass

class UnknownMethodError(TreeTimeError):
    """MissingDataError class raised when tree or alignment are missing"""
    pass

class NotReadyError(TreeTimeError):
    """NotReadyError class raised when results are requested before inference"""
    pass


from .treeanc import TreeAnc
from .treetime import TreeTime, plot_vs_years
from .clock_tree import ClockTree
from .treetime import ttconf as treetime_conf
from .gtr import GTR
from .gtr_site_specific import GTR_site_specific
from .merger_models import Coalescent
from .treeregression import TreeRegression
from .argument_parser import make_parser


