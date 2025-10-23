******************************
Coalescent class documentation
******************************

.. Note::
    Using the coalescent model is optional. When running via the command line this class will only be initialized when
    the flag

    .. code-block:: bash

        --coalescent <arg>

    is used.
    The argument is either ``‘const’`` (have TreeTime estimate a constant coalescence rate), or ``‘skyline’`` (estimate a piece-wise linear merger rate trajectory)
    or a floating point number giving the time scale of coalescence in units of divergence (``Tc``). This is also
    called the effective population size.

.. autoclass::  treetime.Coalescent
    :members: __init__

Parameters of the Kingsman coalescence model
--------------------------------------------

.. automethod:: treetime.Coalescent.branch_merger_rate
.. automethod:: treetime.Coalescent.total_merger_rate
.. automethod:: treetime.Coalescent.cost

Skyline Methods
---------------------

.. automethod:: treetime.Coalescent.optimize_skyline
.. automethod:: treetime.Coalescent.skyline_empirical
.. automethod:: treetime.Coalescent.skyline_inferred