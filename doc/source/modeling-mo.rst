.. _multiple-objectives:


Multiple objectives
----------------------------------

.. image:: images/berg-tal.svg
  :width: 200
  :align: right

To consider multiple objectives in an AMPL model, use
:ref:`solver option <solver-options>` ``obj:multi``.
Otherwise, only the 1st objective is considered
(or any objective specified by ``obj:no``.)

.. code-block:: ampl

    minimize total_cost {s in STORE}:
       sum {j in FOOD} cost[s,j] * Buy[j];

    minimize total_number:  sum {j in FOOD} Buy[j];


See the
`Multi-objective AMPL Colab notebooks <https://colab.ampl.com/tags/multiple-objectives.html>`_
for examples.


Blended objectives
************************************************

By default, all objectives are blended together
(summed up with the corresponding signs.)
Suffixes ``.objweight`` can be used to change the individual weights
and objective senses, according to the option ``obj:multi:weight``.


Lexicographical objectives
********************************************************

To apply hierarchical optimization, use suffix ``.objpriority``,
as described in the ``obj:multi`` option description.

.. code-block:: ampl

    maximize ReverseSeniority {e in 1..2, i in I: E[i]==e}:
      sum {t in V[i]: Pr[i, t]==0}
        S[i] * x[i, t]
      suffix objpriority (2-e)*S_range + 1 + S[i] - min {j in I} S[j];

Suffixes ``.objabstol`` and ``.objreltol`` allow for objective degradation.
