.. _modeling-overview:

Modeling overview and implemented solvers
------------------------------------------


AMPL's newly extended C++ solver interface library, MP, is publicly
available in the `ampl/mp <https://github.com/ampl/mp>`_ repository.
Solver interfaces built with MP are able to handle a significantly
expanded range of model expressions.
Currently available MP-based solvers include:

- `gurobi <https://github.com/ampl/mp/tree/develop/solvers/gurobi>`_,
  an interface to the `Gurobi solver <https://ampl.com/products/solvers/solvers-we-sell/gurobi/>`_

- `cplex <https://github.com/ampl/mp/tree/develop/solvers/cplex>`_,
  an interface to the `IBM ILOG CPLEX solver <https://ampl.com/products/solvers/solvers-we-sell/cplex/>`_

- `copt <https://github.com/ampl/mp/tree/develop/solvers/copt>`_,
  an interface to `Cardinal Optimizer <https://ampl.com/products/solvers/solvers-we-sell/copt/>`_

- `xpress <https://github.com/ampl/mp/tree/develop/solvers/xpress>`_,
  an interface to `FICO Xpress <https://ampl.com/products/solvers/solvers-we-sell/xpress/>`_

- `mosek <https://github.com/ampl/mp/tree/develop/solvers/mosek>`_,
  an interface to the `MOSEK solver <https://ampl.com/products/solvers/solvers-we-sell/mosek/>`_

- `highs <https://github.com/ampl/mp/tree/develop/solvers/highsmp>`_,
  an interface to the open-source `HiGHS solver <https://ampl.com/products/solvers/open-source-solvers/>`_

- `cbc <https://github.com/ampl/mp/tree/develop/solvers/cbcmp>`_,
  an enhanced interface to the open-source
  `CBC solver <https://ampl.com/products/solvers/open-source-solvers/>`_

- `scip <https://github.com/ampl/mp/tree/develop/solvers/scipmp>`_,
  an interface to the open-source `SCIP solver <https://ampl.com/products/solvers/open-source-solvers/>`_

- `gcg <https://github.com/ampl/mp/tree/develop/solvers/gcgmp>`_,
  an interface to the open-source `GCG solver <https://ampl.com/products/solvers/open-source-solvers/>`_

Binaries for these solvers can be downloaded, in distribution
bundles and individually, through the `AMPL Portal <https://portal.ampl.com>`_.
Solver options and features are described
at :ref:`features-guide`
and, for concrete solvers,
at `AMPL Development <https://dev.ampl.com/solvers/index.html>`_.



The expanded MP solver interface library offers new support
for the following categories of operators and expressions:

- Conditional operators: ``if-then-else``; ``==>``, ``<==``, ``<==>``
- Logical operators: ``or``, ``and``, ``not``; ``exists``, ``forall``
- Piecewise linear functions: ``abs``; ``min``, ``max``; ``<<breakpoints; slopes>>``
- Counting operators: ``count``; ``atmost``, ``atleast``, ``exactly``; ``numberof``
- Relational and comparison operators: ``>(=)``, ``<(=)``, ``(!)=``; ``alldiff``
- Complementarity operator: ``complements``
- Nonlinear operators and functions: ``*``, ``/``, ``^``; ``exp``, ``log``;
  ``sin``, ``cos``, ``tan``; ``sinh``, ``cosh``, ``tanh``
- Set membership operator: ``in``

Modeling details and examples are given in the :ref:`expressions_supported` section below.
Technical details, configuration settings, and tools are in the :ref:`modeling-tools` section.
See also the individual solvers' documentation at
`AMPL Development <https://dev.ampl.com/solvers/index.html>`_
for more details of solver-specific features:

- Choice between linearization in the interface and native solver support for some operations
- Handling of AMPL suffixes on constraints that are transformed by the interface

The slides from our presentation on
`Advances in Model-Based Optimization <https://ampl.com/MEETINGS/TALKS/2022_07_Bethlehem_Fourer.pdf>`_
provide overview of the MP interface library in the context of AMPL applications,
including comments on implementation and efficiency issues.

