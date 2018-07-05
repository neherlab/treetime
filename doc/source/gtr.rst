***********************
GTR class documentation
***********************


.. autoclass::  treetime.GTR
    :members: __init__

.. automethod:: treetime.GTR.standard

.. automethod:: treetime.GTR.custom

.. automethod:: treetime.GTR.random

.. automethod:: treetime.GTR.infer

.. automethod:: treetime.GTR.assign_rates

.. Note::
    GTR object can be modified in-place by calling :py:func:`treetime.GTR.assign_rates`


Sequence manipulation
---------------------

.. automethod:: treetime.GTR.compress_sequence_pair

Distance and probability computations
-------------------------------------

.. automethod:: treetime.GTR.optimal_t

.. automethod:: treetime.GTR.optimal_t_compressed

.. automethod:: treetime.GTR.prob_t

.. automethod:: treetime.GTR.prob_t_compressed

.. automethod:: treetime.GTR.prob_t_profiles

.. automethod:: treetime.GTR.propagate_profile

.. automethod:: treetime.GTR.sequence_logLH

.. automethod:: treetime.GTR.expQt
