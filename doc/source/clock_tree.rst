*****************************
ClockTree class documentation
*****************************
ClockTree is a class that implements the core algorithms for maximum likelihood time tree inference. It operates on a tree with fixed topology. All operations the reroot or change tree topology are part of the TreeTime class.

.. .. autoclass:: treetime.ClockTree
..     :members:


ClockTree docstring and constructor
===================================

.. autoclass:: treetime.ClockTree
    :members: __init__


Running TreeTime analysis
=========================

.. automethod:: treetime.ClockTree.init_date_constraints

.. automethod:: treetime.ClockTree.make_time_tree


Post-processing
===============

.. automethod:: treetime.ClockTree.branch_length_to_years

.. automethod:: treetime.ClockTree.convert_dates

.. automethod:: treetime.ClockTree.get_confidence_interval

.. automethod:: treetime.ClockTree.get_max_posterior_region

.. automethod:: treetime.ClockTree.timetree_likelihood
