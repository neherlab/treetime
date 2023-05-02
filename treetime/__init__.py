version="0.10.0"
## Here we define an error class for TreeTime errors, MissingData, UnknownMethod and NotReady errors
## are all due to incorrect calling of TreeTime functions or input data that does not fit our base assumptions.
## Errors marked as TreeTimeUnknownErrors might be due to data not fulfilling base assumptions or due
## to bugs in TreeTime. Please report them to the developers if they persist.
class TreeTimeError(Exception):
    """
    TreeTimeError class
    Parent class for more specific errors
    Raised when treetime is used incorrectly in contrast with `TreeTimeUnknownError`
    `TreeTimeUnknownError` is raised when the reason of the error is unknown, could indicate bug
    """
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

class TreeTimeUnknownError(Exception):
    """TreeTimeUnknownError class raised when TreeTime fails during inference due to an unknown reason. This might be due to data not fulfilling base assumptions or due  to bugs in TreeTime. Please report them to the developers if they persist."""
    pass

import os, sys
recursion_limit = os.environ.get("TREETIME_RECURSION_LIMIT")
if recursion_limit:
    sys.setrecursionlimit(int(recursion_limit))
else:
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))

from .treeanc import TreeAnc
from .treetime import TreeTime, plot_vs_years
from .clock_tree import ClockTree
from .treetime import ttconf as treetime_conf
from .gtr import GTR
from .gtr_site_specific import GTR_site_specific
from .merger_models import Coalescent
from .treeregression import TreeRegression
from .argument_parser import make_parser


