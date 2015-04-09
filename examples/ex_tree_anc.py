 """
Script shows basic functionality of reading the tree and performing basic preparatory operations.

**NOTE** this function shows the  operations of the TreeAnc class, which only purpose is to prepare the tree for the ML date inferrence. The latter is done by the TreeTime class. If you want to see the functionality of this latter class, please refer to the function read_t_tree.
Args:
 - tinf (str): path to tree input file
 - ainf (str): path to alignment input file

"""

import numpy as np
from Bio import AlignIO
from time_tree import tree_anc