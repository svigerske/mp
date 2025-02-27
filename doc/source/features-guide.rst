.. _features-guide:

Features guide for MP-based AMPL solvers
****************************************

.. highlight:: ampl

The MP framework defines standard *solver features* that solvers might support;
these are usually characterized by a set of :ref:`solver-options` used to control the feature,
sometimes suffixes to pass required data and results, and may change the behaviour
of the solution process. Furthermore, MP offers unified :ref:`solve-result-codes`.

This page presents the semantics of the most common solver features; for a development
reference see :ref:`howto`.


.. _solver-options:

Solver options
=================

Solver options are key-value pairs controlling a solver's behavior.
We distinguish between :ref:`ampl-solver-options` and :ref:`native-options`.


.. _ampl-solver-options:

AMPL/MP solver options
-----------------------------

AMPL/MP solver options provide a unified interface to
MP parameters, as well as underlying solver's configuration:

.. code-block:: ampl

    ampl: option solver gurobi;
    ampl: option gurobi_options 'outlev=1';      ## verbose output
    ampl: solve;

Many of them are
`standardized across AMPL solvers <https://dev.ampl.com/solvers/index.html>`_.


List all available options
'''''''''''''''''''''''''''''''''''''''''''''

Run the AMPL solver executable with ``-=`` to list all options,
or with ``-=key`` to list all options containing ``key``:

.. code-block:: bash

    $ highs -=acc:
    .....
    acc: Options:

    acc:linrange (acc:linrng)
        Solver acceptance level for 'LinConRange', default 2

On the Web, options for AMPL solvers are published in 
`AMPL Development <https://dev.ampl.com/solvers/index.html>`_.


Option file
'''''''''''''''''''''''''''''''''''

It is possible to input a file with predefined AMPL solver options,
using ``tech:optionfile``:

.. code-block:: ampl

    ampl: option cbc_options 'optionfile="options_experiment1.txt"';

For the underlying solver's native options see :ref:`native-options`.


Solver-specific vs common MP options
''''''''''''''''''''''''''''''''''''''''''''''''

From AMPL, options can be passed to an MP solver in two ways.
Solver-specific options are passed
via AMPL option ``(solvername)_options``:

.. code-block:: ampl

    ampl: option solver gurobi;
    ampl: option gurobi_options 'iis=1';         ## Tell Gurobi to find IIS
    ampl: solve;
    Gurobi 10.0.2:   alg:iisfind = 1
    Gurobi 10.0.2: infeasible problem

If your solver executable has a different name, say ``gurobi25.exe``,
and you provide ``gurobi25_options``, they will be used instead.

Common options (for all MP solvers) can be passed via AMPL option
``mp_options``:

.. code-block:: ampl

    ampl: option mp_options 'lim:time=300';      ## Options for all MP solvers

The value of ``mp_options`` is parsed before ``(solvername)_options``.
Thus, ``mp_options`` allows setting parameters for all MP solvers,
with the possibility to override some of them for a specific solver.



Set options from command line
'''''''''''''''''''''''''''''''''''''''''''''''

When running from command line, there are two ways to pass options:
via the environment variable, or via arguments:

.. code-block:: bash

    mp_options='outlev=1 tech:writeprob=model.mps' gurobi.exe model.nl  ## Method 1
    gurobi.exe model.nl outlev=1 tech:writesol=model.sol                ## Method 2


Query option values
''''''''''''''''''''''''''''''''''''''

To query the value of an option (default, or set via other methods),
use '?' as argument:

.. code-block:: ampl

    ampl: option mosek_options 'threads=?';
    ampl: solve;
    MOSEK 10.0.43:   tech:threads = 0
    MOSEK 10.0.43: optimal; objective 10



.. _native-options:

Native solver options
-------------------------

Some AMPL solvers allow passing native options to the underlying solver.
For example, to control Gurobi ``NumericFocus`` setting, there are two ways:

.. code-block:: ampl

    ampl: option gurobi_options 'alg:numericfocus 3';             ## standard way
    ampl: option gurobi_options 'tech:optionnative "NumericFocus 3"';   ## native
    gurobi.exe model.nl optnative="numericfocus 2" optnative="Seed 500" # cmdline

Additionally, for some solvers, native options can be read / written
from / to files using ``tech:optionnativeread`` and ``tech:optionnativewrite``.


.. _solve-result-codes:

Solve result codes
=================================

The result of the last solve in AMPL can be seen as follows.

.. code-block:: ampl

    ampl: model party2.mod
    ampl: data party2.dat
    ampl: display solve_result_num, solve_result;
    solve_result_num = -1
    solve_result = '?'

    ampl: option solver gurobi;
    ampl: option gurobi_options 'lim:time 20';
    ampl: solve;
    Gurobi 11.0.0:   lim:time = 20
    Gurobi 11.0.0: time limit, feasible solution
    35671 simplex iteration(s)
    1 branching node(s)
    absmipgap=27, relmipgap=0.818182
    ampl: display solve_result_num, solve_result;
    solve_result_num = 402
    solve_result = limit

    ampl: option solve_result_table;
    option solve_result_table '\
    0       solved\
    100     solved?\
    200     infeasible\
    300     unbounded\
    400     limit\
    500     failure\
    ';

MP details the solve result codes as follows:

.. code-block:: ampl

    ampl: shell "mosek -!";
    Solve result table for MOSEK 10.2.0
        0- 99	solved: optimal for an optimization problem,
              feasible for a satisfaction problem
      100-199	solved? solution candidate returned but error likely
          150	solved? MP solution check failed (option sol:chk:fail)
      200-299	infeasible
      300-349	unbounded, feasible solution returned
      350-399	unbounded, no feasible solution returned
      400-449	limit, feasible: stopped, e.g., on iterations or Ctrl-C
      450-469	limit, problem is either infeasible or unbounded.
              Disable dual reductions or run IIS finder for definitive answer.
      470-499	limit, no solution returned
      500-999	failure, no solution returned
          550	failure: numeric issue, no feasible solution

Individual solvers may add more specific values in the corresponding ranges.
To list solver-specific codes, use command-line switch ``-!`` as above,
or visit `AMPL Development <https://dev.ampl.com/solvers/index.html>`_.
More information is in Chapter 14 of the
`AMPL Book <https://ampl.com/learn/ampl-book/>`_.
See also the roll cutting example on `AMPL Colab <https://colab.ampl.com>`_.

Solvers support
===============

.. |y| unicode:: U+2705 
   :trim:
  
.. |n| unicode:: U+274C
   :trim:

.. |e| unicode:: U+2713
   :trim:

This table summarizes the solver driver support for each solver feature; some features, 
denoted by |e| are supported through emulation (aka they were not available in the original 
solver but are emulated by MP). 

.. _Copt: https://dev.ampl.com/solvers/copt/
.. _CPLEX: https://dev.ampl.com/solvers/cplex/
.. _Gurobi: https://dev.ampl.com/solvers/gurobi/
.. _Mosek: https://dev.ampl.com/solvers/mosek/
.. _Xpress: https://dev.ampl.com/solvers/xpress/
.. _CBC: https://dev.ampl.com/solvers/cbc/
.. _GCG: https://dev.ampl.com/solvers/gcg/
.. _HiGHS: https://dev.ampl.com/solvers/highs/
.. _SCIP: https://dev.ampl.com/solvers/scip/

+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| Feature                      | Copt_ | CPLEX_ | Gurobi_ | Mosek_ | Xpress_ | CBC_ | GCG_ | HiGHS_ | SCIP_ |
+==============================+=======+========+=========+========+=========+======+======+========+=======+
| :ref:`writeprob`             | |y|   | |y|    | |y|     | |y|    | |y|     | |y|  | |y|  | |y|    | |y|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| Solution Export              | |n|   | |y|    | |y|     | |n|    | |y|     | |n|  | |n|  | |y|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`report-times`          | |y|   | |y|    | |y|     | |y|    | |y|     | |y|  | |y|  | |y|    | |y|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`feasibilityrelaxation` | |y|   | |y|    | |y|     | |n|    | |n|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`iis`                   | |y|   | |y|    | |y|     | |n|    | |y|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`basisio`               | |y|   | |y|    | |y|     | |y|    | |y|     | |y|  | |y|  | |y|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`kappa`                 | |n|   | |y|    | |y|     | |n|    | |n|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`multiplesolutions`     | |y|   | |y|    | |y|     | |n|    | |y|     | |n|  | |y|  | |n|    | |y|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`multipleObjectives`    | |e|   | |y|    | |y|     | |e|    | |y|     | |e|  | |e|  | |e|    | |e|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`sensitivityAnalysis`   | |n|   | |n|    | |y|     | |y|    | |n|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`unboundedRays`         | |y|   | |y|    | |y|     | |y|    | |n|     | |n|  | |n|  | |y|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`warm-start`            | |y|   | |y|    | |y|     | |y|    | |y|     | |y|  | |n|  | |y|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`fixedModel`            | |n|   | |y|    | |y|     | |n|    | |y|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`lazyConstraints`       | |n|   | |n|    | |y|     | |n|    | |n|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`returnMIPgap`          | |y|   | |y|    | |y|     | |y|    | |y|     | |n|  | |y|  | |y|    | |y|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`returnBestBound`       | |y|   | |y|    | |y|     | |y|    | |y|     | |n|  | |y|  | |y|    | |y|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+
| :ref:`varPriorities`         | |n|   | |n|    | |y|     | |n|    | |n|     | |n|  | |n|  | |n|    | |n|   |
+------------------------------+-------+--------+---------+--------+---------+------+------+--------+-------+


* |n| = Not supported
* |e| = Emulated
* |y| = Supported natively



General features
================


.. _outlev:

Output level
------------

All solvers print some information during the solution process. For the sake of clarity,
when using the solvers from AMPL the information is mostly ommitted, and only the ``solver message`` is
displayed, which describes the outcome of the optimization process.
The amount of information printed is controlled by the option ``outlev``::

    option <solver>_options 'outlev=1';

In case of repeated executions it is often desiderable to hide all solver output. The solver message 
can be suppressed by setting the option ``solver_msg`` to 0. The solver banner is displayed anyway; 
to suppress it a redirection is needed::

  option solver_msg 0;
  solve > NUL;

*Related option*: `tech:timing`.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``outlev``
   * - **Applicability**
     - All models
   * - **Values**
     - Values:

       * **0** - No (default)
       * **1** - Yes, print detailed execution log
       * **2, 3, ...** - Solver-dependent values


.. _writeprob:

Model export
------------

Most solvers can export the model before solving. This is usually
controlled by the option ``writeprob``::

    option <solver>_options 'writeprob=/tmp/diet.lp';


The format is solver-dependent and determined by the file extension
('.lp' in the example).


.. list-table::
   :header-rows: 0

   * - **Option**
     - ``writeprob``
   * - **Applicability**
     - All models
   * - **Values**
     - Values:

       * **filename** - Filename for the exported model


.. _report-times:

Report solution time
--------------------

All solvers can measure the execution times of the solution process;
is this option is sset to 1 the solver will display and return via
problem suffixes this information. Setting it to 2 will give more granular
information. The reported times with the options set to 1 are: `time_solver` 
(the time taken by the solver itself to solve the model), `time_setup` (the time spent in the solver driver for
reading the model file, reformulating it and prepare all the memory structures)
and  `time` (=time_solver+time_setup, the total time spent in the driver). By setting
it to 2 `time_read`, `time_conversion` and `time_output` are also reported. Example::


    option <solver>_options 'timing=1';


.. list-table::
   :header-rows: 0

   * - **Option**
     - ``tech::timing``, ``tech:reporttimes``
   * - **Applicability**
     - All models
   * - **Input**
     - None
   * - **Output**
     - Suffixes (on problem):

       * ``time`` total time spent in the driver
       * ``time_setup`` time spent to read, reformulate and allocate the model
       * ``time_solver`` time spent solving the model
       * ``time_read`` time spent to read the NL file
       * ``time_conversion`` time spent in the reformulation process
       * ``time_output`` time spent reporting the results
   
   * - **Values**
     - Values:
  
       * **0** - No (default)
       * **1** - Yes, basic information
       * **2** - Yes, granular information 

   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         option <solver>_options "tech:reporttimes=1";
         solve;

         display Initial.time, Initial.time_setup, Initial.time_solver;

       Output:

       .. code-block:: bash

          tech:reporttimes = 1
          Setup time:    0.0113952s
          Solution time: 0.0016453s

          suffix time OUT;
          suffix time_setup OUT;
          suffix time_solver OUT;

          Initial.time = 19.43
          Initial.time_setup = 0.88
          Initial.time_solver = 19.46


.. _warm-start:

Warm start
----------

Solution process can often benefit of a solution (a set of variable values) to start the algorithm. 
This is passed to supporting solver automatically if the option is activated and variables in AMPL
have a value assigned. Note that, for LP problems, also the dual values can be passed.

This option controls whether to use incoming primal (and dual, for LP) variable values in 
a warmstart.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:start``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - Variable values
   * - **Output**
     - None
   * - **Values**
     - Sum of:

       * **0** - No (default)
       * **1** - Yes (for LP: if there is no incoming alg:basis) (default)
       * **2** - Yes (for LP: ignoring the incoming alg:basis, if any)
   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         option <solver>_options "alg:start=1";
         let n:=250; # increase the size of the model to have noticeable solution times

         solve;
         printf "Solution without warm start took %fs\n", _solve_time;

         # Now an optimal solution is already present, we pass it to the solver
         # which would use to start the solution process
         solve;
         printf "Solution with warm start took %fs\n", _solve_time;

       Output:

       .. code-block:: shell
 
            ...

            Solution without warm start took 2.89062s
            
            ...

            Solution with warm start took 0.671875s


.. _basisio:

Input and output basis
----------------------

A basis is a set of variable values representing a feasible and extreme solution.
Simplex solvers normally calculate this as part of the solution process, while
interior point methods must perform additional steps (crossover) to get it.
In a way similar to :ref:`warm start <warm-start>`, a basis can also be passed to the solver,
which will use it as starting point for searching for a solution.

This option controls whether to use or return a basis.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:basis``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - Suffixes:

       * ``sstatus`` on variables and constraints
   * - **Output**
     - Suffixes:

       * ``sstatus`` on variables and constraints
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Use incoming basis (if provided)
       * **2** - Return final basis
       * **3** - Both (1 + 2, default)

   * - **Example**
     - Use :ref:`this model <multiObjectiveDiet>`

       Execute::

         option gurobi_options "alg:start=0 outlev=1"; # disable passing the solution

         solve;
         display Buy.sstatus; # display basis status

         solve; # second solve with take much less although a solution is not provided

       In the solver logs, we can see the expected behaviour:

       .. code-block:: shell

          x-Gurobi 9.5.2: optimal solution; objective 74.27382022
          3 simplex iterations
          Objective = total_cost['A&P']
          ampl: display Buy.sstatus;
          Buy.sstatus [*] :=
          BEEF  low
          CHK  upp
          FISH  low
          HAM  low
          MCH  low
          MTL  bas
          SPG  bas
          TUR  low;

          ... # second solve:          

          Solved in 0 iterations and 0.00 seconds (0.00 work units)
          

.. _feasibilityrelaxation:

Feasibility Relaxation
----------------------

The feasibility relaxation functionality enables the solver to find a feasible
solution even if the original model is unfeasible without explicitly adding
slack variables to the constraints.
In the feasibility relaxation problem, 

#. Each variable :math:`x` can violate its bounds (:math:`lb \leq x \leq ub`):
  
   * Violation of lower bound :math:`lbv = max(0, lb-x)`
   * Violation of upper bund :math:`ubv = max(0, x-ub)`

#. Each constraint body :math:`c` can violate its bounds also (:math:`c \leq rhs`)

   * Constraint violation :math:`rhsv = max(0, c-rhs)`

The objective then becomes to minimize some function of the
violations (e.g. the number of violations, or their sum - possibly weighted by some
penalty values).
The penalty values (used in some kinds of feasibility relaxation problmes) can be 
controlled with macro defaults (e.g. option ``alg:ubpen`` sets the penalty weight for 
all upper buonds violations, and its default values is 1) or, with more granularity,
on each entity via suffix values (e.g. variable suffix ``ubpen`` on variables, default
value 0). 
Penaly weights < 0 are treated as Infinity, allowing no violation.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:feasrelax``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     -  * Options
  
          * ``alg:lbpen``: penalty for lower bound violations if suffix ``lbpen`` is not defined - default 1
          * ``alg:ubpen``: penalty for upper bound violations if suffix ``ubpen`` is not defined - default 1
          * ``alg:rhspen``: penalty for rhs violations if suffix ``rhspen`` is not defined - default 1

        * Suffixes

          * ``lbpen`` on variables - penalty for lower bound violations - default 0
          * ``ubpen`` on variables - penalty for upper bound violations - default 0
          * ``rhspen`` on constraints - penalty for rhs violations - default 0
   * - **Output**
     - None
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Yes, minimizing the weighted sum of violations
       * **2** - Yes, minimizing the weighted sum of squared violations
       * **3** - Yes, minimizing the weighted count of violations
       * **3-6** - Same objective as 1-3, but also optimize the original objective, subject to the violation being minimized
   * - **Example**
     - Use :ref:`this model <infeasibleModel>`

       Solve the model changing the penalties to get different solutions::

          options <solver>_options "alg:feasrelax=1";
          option presolve 0;    # otherwise the model could be oversimplified

          solve; display x,y;

       Gives:

        .. code-block:: bash

          x = 1
          y = 1
        
       Now we want to force all variable lower bounds to be respected::

          options <solver>_options "alg:feasrelax=1 alg:lbpen=-1";
          solve; display x,y;

       Gives, as expected x=5 (it had a lower bound of 5):

        .. code-block:: bash

            x = 5
            y = 1
        
       Single violations can be controlled via suffixes; in this case we
       want to control the constraints violations::

          options <solver>_options "alg:feasrelax=1 alg:lbpen=1"; # allow lower bounds to be violated again

          suffix rhspen IN;
          let C1.rhspen := 1; # normal weight
          let C2.rhspen := -1; # C2 can NOT be violated
          let C3.rhspen := 10; # We'd rather not violate C3
          solve;
          display C1.slack, C2.slack, C3.slack;

       Gives:

        .. code-block:: bash

          C1.slack = -3
          C2.slack = 4
          C3.slack = 0

       C1 is violated - which makes sense as we specified that C2 cannot be violated
       and we gave an higher avoidance weight to C2. If we want to violate C1 instead,
       we can::

          let C1.rhspen := 10; # We'd rather not violate C1
          let C2.rhspen := 1; # C2 can be violated
          let C3.rhspen := 1; # Normal weight
          solve;
          display C1.slack, C2.slack, C3.slack;

       Which gives, as expected:

        .. code-block:: bash

          C1.slack = 18
          C2.slack = 2
          C3.slack = -1


.. _multiplesolutions:

Multiple solutions
------------------

More often than not, optimization problems have more than one optimal solution; moreover, during the 
solution process, MIP solvers usually find sub-optimal solutions, which are normally discarded.
They can be however be kept, and in most cases there are solver-specific options to control how
the search for additional solutions is performed.

The main (and generic) options that controls the search are ``sol:stub`` amd ``sol:count``, which
control respecitvely the base-name for the files where additional solution will be stored and
if to count additional solutions and return them in the ``nsol`` problem suffix.
Specifying a stub name automatically enables the solutions count; found solutions are written to 
files [``solutionstub1.sol'``,  ... ``solutionstub<nsol>.sol``].


.. list-table::
   :header-rows: 0

   * - **Option**
     - ``sol:stub``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffixes ``nsol`` and ``npool`` on problem
   * - **Values**
     - The name used as base file name for the alternative solutions

   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         option gurobi_options "sol:stub=queentake";
         solve;
         printf "I have found %d solutions\n", Initial.nsol;

         printf "Displaying solution 1";
         solution queentake1.sol;
         display X;

         printf "Displaying solution 2";
         solution queentake2.sol;
         display X;

       Output:

       .. code-block:: shell
 
            x-Gurobi 9.5.2: sol:stub=queentake
            x-Gurobi 9.5.2: optimal solution; objective 10
            355 simplex iterations
            23 branching nodes

            suffix nsol OUT;
            suffix npool OUT;

            I have found 10 solutions
            Displaying solution 1
            x-Gurobi 9.5.2: Alternative solution; objective 10
            X [*,*]
            :    1   2   3   4   5   6   7   8   9  10    :=
            1    0   0   0   0   0   0   1   0   0   0
            2    1   0   0   0   0   0   0   0   0   0
            3    0   0   0   1   0   0   0   0   0   0
            4    0   0   0   0   0   1   0   0   0   0
            5    0   0   0   0   0   0   0   0   1   0
            6    0   0   1   0   0   0   0   0   0   0
            7    0   0   0   0   0   0   0   0   0   1
            8    0   0   0   0   0   0   0   1   0   0
            9    0   1   0   0   0   0   0   0   0   0
            10   0   0   0   0   1   0   0   0   0   0;

            Displaying solution 2
            x-Gurobi 9.5.2: Alternative solution; objective 10
            X [*,*]
            :    1   2   3   4   5   6   7   8   9  10    :=
            1    0   0   1   0   0   0   0   0   0   0
            2    0   0   0   0   0   1   0   0   0   0
            3    0   0   0   0   0   0   0   0   1   0
            4    1   0   0   0   0   0   0   0   0   0
            5    0   0   0   1   0   0   0   0   0   0
            6    0   0   0   0   0   0   1   0   0   0
            7    0   0   0   0   0   0   0   0   0   1
            8    0   0   0   0   0   0   0   1   0   0
            9    0   1   0   0   0   0   0   0   0   0
            10   0   0   0   0   1   0   0   0   0   0


.. _sensitivityAnalysis:

Sensitivity analysis
--------------------

It is often useful to know the ranges of variables and constraint bodies for which the optimal basis
remains optimal. Solvers supporting this feature return such ranges in suffixes after solving to optimum.
This option controls whether to calculate these values and return them in the suffixes listed below.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:sens``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffixes for variables only:

       * ``sensobjlo`` smallest objective coefficient
       * ``sensobjhi`` greatest objective coefficient

       Suffixes for variables and constraints:

       * ``senslblo`` smallest variable/constraint lower bound
       * ``senslbhi`` greatest  lower bound
       * ``sensublo`` smallest  upper bound
       * ``sensubhi`` greatest  upper bound

   * - **Values**
     - Sum of:

       * **0** - No (default)
       * **1** - Yes
       
   * - **Example**
     - Use :ref:`multiObjectiveDiet`

       Execute::

          options <solver>_options "alg:sens=1"; 
          option presolve 0;  # disable model transformations in AMPL
          solve;

       Then the ranges for variables and constraints can be examined::

          display Buy, Buy.sstatus, Buy.sensublo, Buy.sensubhi;

       Which gives:

       .. code-block:: bash

          display Buy.sensublo, Buy.sensubhi, Buy.sstatus, Buy;
          :    Buy.sensublo  Buy.sensubhi Buy.sstatus     Buy       :=
          BEEF    2           1e+100        low          2
          CHK     9.12987         10.9792   upp         10
          FISH    2           1e+100        low          2
          HAM     2           1e+100        low          2
          MCH     2           1e+100        low          2
          MTL     6.23596     1e+100        bas          6.23596
          SPG     5.25843     1e+100        bas          5.25843
          TUR     2           1e+100        low          2;          

.. _kappa:

Kappa
-----

Kappa is the condition number for the current LP basis matrix. 
It is a measure of the stability of the current solution :math:`Ax=b`
measuring the rate at which the solution :math:`x` will change with respect to a 
change in :math:`b`. 
It is only available for basic solutions, therefore it is not available for barrier method
if crossover is not applied.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:kappa``
   * - **Applicability**
     - LP and MIP models with optimal basis
   * - **Input**
     - None
   * - **Output**
     - Additional text in ``solve_message`` and suffix:

       * ``kappa`` on objective and problem
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Report kappa in solve_message
       * **2** - Return kappa in the solver-defined suffix ``kappa``
   * - **Example**
     - Use :ref:`multiObjectiveDiet`

       Solve the model and report kappa::

          options <solver>_options "alg:kappa=3"; 
          solve;

          display Initial.kappa, total_number.kappa;
 
       Gives:

       .. code-block:: bash

        x-Gurobi 9.5.2: optimal solution; objective 30.92537313
        kappa value: 53.5399
        5 simplex iterations

        suffix kappa OUT;

        Initial.kappa = 53.5399
        total_number.kappa = 53.5399


.. _unboundedRays:

Unbounded rays
--------------

When a model is unbounded, a vector :math:`r` (unbounded ray)  can be found such that
when added to any feasible solution :math:`x`, the resulting vector is a feasible solution
with an improved objective value.
When a model is infeasible, the dual solution is unbounded and the same as above can be applied
to constraints.
This option controls whether to return suffix ``unbdd`` if the objective is unbounded
or suffix ``dunbdd`` if the constraints are infeasible.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:rays``
   * - **Applicability**
     - Unbounded/unfeasible LP and MIP linear models
   * - **Input**
     - None
   * - **Output**
     - Suffixes:

       * ``unbdd`` on variables if the problem is unbounded
       * ``dunbdd`` on constraints if the problem is infeasible
   * - **Values**
     - Sum of:

       * **0** - Do not calculate or return unbounded rays
       * **1** - Return only ``unbdd``
       * **2** - Return only ``dunbdd``
       * **3** - Return both (default)
   * - **Example**
     - Use :ref:`this model <iisModel>`

       Solve the (infeasible) model and report the ``dunbdd`` rays::

          options gurobi_options "alg:rays=3"; # it is default already
          option presolve 0;  # else the model could be oversimplified
          solve;
          display c1.dunbdd, c2.dunbdd, c3.dunbdd;

       Output:

       .. code-block:: shell

          x-Gurobi 9.5.2: alg:rays=3
          x-Gurobi 9.5.2: infeasible problem

          suffix dunbdd OUT;

          c1.dunbdd = 1
          c2.dunbdd = 0
          c3.dunbdd = 0


.. _multipleObjectives:

Multiple objectives
-------------------

Many real world problems have multiple objectives; often this scenario is tackled by blending all the objectives
by linear combination when formulating the model, or by minimizing each unwanted objective deviations from a pre-specified
goal.
Many solvers can facilitate the formulation; the available functionalities are solver-specific. For other solvers,
MP :ref:`emulates the multi-objective capability <multiple-objectives>`. Consult the ``obj:multi``
:ref:`option <solver-options>` documentation
for the functionalities available on your solver.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``obj:multi``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - None
   * - **Output**
     - All objectives are reported in ``solve_message``

   * - **Values**
     - Values:

       * **0** - No (default)
       * **1** - Yes, natively supported if available, otherwise emulated
       * **2** - Yes, emulated
   * - **Example**
     - Use :ref:`multiObjectiveDiet`

       Execute::

          options <solver>_options "obj:multi=1"; 
          solve;

       Output:

       .. code-block:: shell

          x-Gurobi 9.5.1: obj:multi=1
          x-Gurobi 9.5.1: optimal solution; objective 74.27382022
          Individual objective values:
            _sobj[1] = 74.27382022
            _sobj[2] = 75.01966292
            _sobj[3] = 79.59719101
            _sobj[4] = 31.49438202


   
Irreducible Inconsistent Subset (IIS)
------------------------------------------

Given an infeasible model, it is useful to know where the infeasibility comes from, that is, which
bounds and/or constraints are incompatible.
An IIS is a subset of constraints and variables that is still infeasible and where if a member is removed
the subsystem becaomes feasible. Note that an infeasible model can have more than one IIS.

This options controls whether to perform the additional computational steps required to find an IIS.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:iisfind``, ``iisfind``,  ``iis``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffix ``iis`` on variables and constraints

   * - **Values**
     - Values:

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Use :ref:`this model <iisModel>`

       Execute::

        option <solver>_options "iisfind=1";
        option presolve 0;  # else the model could be oversimplified
        solve;
        display c1.iis, c2.iis, c3.iis;

       Output:

       .. code-block:: shell

          x-Gurobi 9.5.2: alg:iisfind=1
          x-Gurobi 9.5.2: infeasible problem
          2 simplex iterations

          suffix iis symbolic OUT;
          suffix dunbdd OUT;
  
          c1.iis = mem
          c2.iis = mem
          c3.iis = mem


MIP-only features
=================


.. _returnMIPgap:

Return MIP gap
--------------

The MIP gap describes how far the reported solution for the MIP model is from the
best bound found. It is defined as:

:math:`absgap = | objective - bestbound |`

:math:`relmipgap = absgap / | objective |``

It gives a measure of how far the current solution is from the
theoretical optimum.
The AMPL option controls whether to return mipgap suffixes or include mipgap values 
in the solve_message. Returned suffix values are ``+Infinity`` if no integer-feasible 
solution has been found, in which case no mipgap values are reported in the solve_message.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:return_gap``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None
   * - **Output**
     - Additional text in ``solve_message`` and suffixes:

       * ``relmipgap`` on objective and problem
       * ``absmipgap`` on objective and problem
   * - **Values**
     - Sum of:

       * **0** - Do not report gaps
       * **1** - Return .relmipgap suffix (relative to ``|obj|``)
       * **2** - Return .absmipgap suffix (absolute mipgap)
       * **4** - Suppress mipgap values in solve_message
   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         # 3 = return both absolute and relative MIP gap, and also report them
         # in the solve_message
         option <solver>_options "return_mipgap=3";
         solve;

         display max_queens.relmipgap, max_queens.absmipgap;
         display Initial.relmipgap, Initial.absmipgap;

       Output:

       .. code-block:: bash

          x-Gurobi 9.5.2: mip:return_gap=3
          x-Gurobi 9.5.2: optimal solution; objective 10

          suffix absmipgap OUT;
  
          max_queens.relmipgap = 0
          max_queens.absmipgap = 0
          Initial.relmipgap = 0
          Initial.absmipgap = 0



.. _returnBestBound:

Return best dual bound
----------------------

The best dual bound (on the objective) represents what is the currently 
proven best value that the objective value can assume. Usually solvers terminate
when the current solution is close enough to the best bound.

This option controls whether to return suffix .bestbound for the best known MIP dual
bound on the objective value. The returned value is -Infinity for minimization
problems and +Infinity for maximization problems if there are no integer 
variables or if a dual bound is not available.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:bestbound``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffix:

       * ``bestbound`` on objective
   * - **Values**
     - Sum of:

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

            option <solver>_options "mip:bestbound=3";
            solve;
            display max_queens.bestbound, Initial.bestbound;

       Output:

       .. code-block:: shell

            x-Gurobi 9.5.2: mip:bestbound=1
            x-Gurobi 9.5.2: optimal solution; objective 10
            
            max_queens.bestbound = 10
            Initial.bestbound = 10

.. _lazyConstraints:

Lazy constraints and user cuts
------------------------------

The solution process of a MIP model can be helped by further specifying its structure.
Specifically, constraints can be marked as ``lazy`` or as ``user cuts``.
Such constraints are not included initially, then:

``Lazy constraints`` are pulled in when a feasible solution is found; if the solution violates
them, it is cut off. They are integral part of the model, as they can cut off integer-feasible 
solutions.

``User cuts`` are pulled in also, but they can only cut off relaxation solutions; they are
consider redundant in terms of specifying integer feasibility.

This option controls whether to recognize the suffx ``lazy`` on constraints, which should then be 
positive to denote a lazy constraint and negative to mark a contraint as a user cut.


.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:lazy``
   * - **Applicability**
     - MIP models
   * - **Input**
     - Suffix:

       * ``lazy`` on constraints (>0 for lazy constraint, <0 for user cuts) 
   * - **Output**
     - None
  
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Accept >0 values to denote lazy constraints
       * **2** - Accept <0 values to denote user cuts
       * **3** - Accept both (default)
   * - **Example**
     - Following ampl model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            option presolve 0;  # disable model transformations in AMPL

            suffix lazy IN;
            let c1.lazy := 1;  # lazy constraint
            let c2.lazy := -1; # user cut, must be redundant as it might never be pulled in
            solve;

            # TODO SHOW OUTPUT



.. _varPriorities:

Variable priorities
-------------------

Solution of MIP models via branch and bound can often be helped by providing
preferences on which variables to branch on. Those can be specified in AMPL via the suffix 
``priority``.

This option controls whether to read those values and use them in the solution process.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:priorities``
   * - **Applicability**
     - MIP models
   * - **Input**
     - Suffix:

       * ``priority`` on variables 
   * - **Output**
     - None
   * - **Values**
     - Values:

       * **0** - Ignore priorities
       * **1** - Read priorities (default)
   * - **Example**
     - Following AMPL model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            option presolve 0;  # disable model transformations in AMPL

            let x.priority := 1;
            let y.priority := 5;
            solve;

            # TODO SHOW OUTPUT




.. _fixedModel:

Fixed model (return basis for MIP)
----------------------------------

At the end of the solution process for a MIP model, a continuous relaxation of the model
with all the integer variables fixed at their integer-optimum value. Some continuous variables
can also be fixed to satisfy SOS or general constraints.
The model can therefore be solved without these types of restrictions to calculate a basis,
dual values or sensitivity information that wouldn't normally be available for MIP problems.

This option controls if to generate and solve the fixed model after solving the integer problem.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:basis``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None

   * - **Output**
     - Suffixes:
  
       * ``dual`` on variables
       * ``sstatus`` on variables
       * See :ref:`sensitivityanalysis` if requested
  
   * - **Values**
     - Values

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Following ampl model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            option <solver>_options "mip:basis=1";

            # TODO SHOW OUTPUT


* Round


Reference models
================
.. toctree::
   :maxdepth: 2

    Reference models <models>

