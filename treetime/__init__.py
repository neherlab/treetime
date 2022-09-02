version="0.9.2"
## Here we define an error class for TreeTime errors, MissingData, UnknownMethod and NotReady errors
## are all due to incorrect calling of TreeTime functions or input data that does not fit our base assumptions.
## Errors marked as TreeTimeOtherErrors might be due to data not fulfilling base assumptions or due 
## to bugs in TreeTime. Please report them to the developers if they persist. 
class TreeTimeError(Exception):
    """TreeTimeError class"""
    pass

class MissingDataError(TreeTimeError):
    """MissingDataError class raised when tree or alignment are missing"""
    pass

class UnknownMethodError(TreeTimeError):
    """MissingDataError class raised when an unknown method is called"""
    pass

class NotReadyError(TreeTimeError):
    """NotReadyError class raised when results are requested before inference"""
    pass

class TreeTimeOtherError(TreeTimeError):
    """TreeTimeOtherError class raised when TreeTime fails during inference due to an unknown reason. This might be due to data not fulfilling base assumptions or due  to bugs in TreeTime. Please report them to the developers if they persist."""
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


