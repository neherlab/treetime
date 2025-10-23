****************************
TreeTime class documentation
****************************
TreeTime is the top-level wrapper class of the time tree inference package. In addition to inferring time trees, TreeTime can reroot your tree, resolve polytomies, mark tips that violate the molecular clock, or infer coalescent models. The core time tree inference is implemented in the class ClockTree.

TreeTime docstring and constructor
==================================

.. autoclass:: treetime.TreeTime
    :members: __init__



Main pipeline method
====================

.. automethod:: treetime.TreeTime.run

Additional functionality
========================

.. automethod:: treetime.TreeTime.resolve_polytomies

.. automethod:: treetime.TreeTime.relaxed_clock

.. automethod:: treetime.TreeTime.clock_filter

.. automethod:: treetime.TreeTime.reroot

.. automethod:: treetime.TreeTime.plot_root_to_tip

.. automethod:: treetime.TreeTime.print_lh

.. autofunction:: treetime.plot_vs_years

