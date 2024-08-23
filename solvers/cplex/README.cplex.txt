The solver "cplex" uses CPLEX (a trademark of IBM,
Inc.; see https://www.ibm.com/products/ilog-cplex-optimization-studio) 
to solve integer, mixed-integer, and linear programming problems; 
it is an alternative to the ASL-based driver implemented using the mp 
library (https://github.com/ampl/mp) for communicating with AMPL and 
for reformlation of certain types of problems.
Normally gurobi is invoked by AMPL's solve command, which gives the 
invocation

     cplex stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
cplex writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver cplex;
     solve;

You can control cplex by setting the environment variable cplex_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $cplex_options.  To see
the possibilities, invoke

        cplex -=

----------
INSTALLING
==========

On Linux systems, libcplex*.so (where the value of "*" depends
on the current version of CPLEX) and the libcplex.so.* to which
it points need to appear in the current directory when CPLEX
itself appears there, or in one of the standard places (specified by
/etc/ld.so.conf on some systems), or in a directory named in
$LD_LIBRARY_PATH.  An alternative is to add a short shell script,
such as

        #!/bin/sh
        LD_LIBRARY_PATH=/usr/local/lib
        export LD_LIBRARY_PATH
        exec /usr/local/bin/cplexx "$@"

to a directory in your usual $PATH (and mark the script executable
with, e.g., "chmod +x cplex").  The above script assumes that the
true "cplex" binary has been moved to /usr/local/bin/cplexx and that
the libcplex* files have been moved to /usr/local/lib.

MacOSX systems are similar to Linux systems, but with DYLD_LIBRARY_PATH
in place of LD_LIBRARY_PATH. 

On MS Windows systems, cplex.exe and the relevant cplex*.dll must
appear somewhere in your usual search $PATH (or in the current
directory).

If you have questions about or find bugs with this stuff,
please contact:

     AMPL Support
     support@ampl.com
