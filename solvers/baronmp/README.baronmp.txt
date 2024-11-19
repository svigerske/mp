Baron Solver
============

Solver baron is a global nonlinear solver based on the BARON solver by
Nikolaos V. Sahinidis and co-authors.  See

        http://www.minlp.com/

For much more detail on BARON, including pointers to papers, see the
accompanying baron_manual.pdf.

It is most convenient to invoke baron with AMPL's solve command:

        ampl: option solver baron;
        ampl: solve;

but baron can also be run separately, with invocation syntax

        baron stub [-AMPL] [keywd=[value] ...]

in which stub comes from AMPL's write command.  For example,

        > ampl -obfoo foo.mod foo.dat
        > baron foo

demonstrates running baron separately; the ampl invocation writes file
foo.nl (stub = "foo"), which baron reads.

Command-line arguments to baron either have the form

        keywd
or
        keywd=value

where keywd is one of the keywords described below.  Alternatively,
you can invoke baron the way AMPL's solve command does, i.e.,

        baron stub -AMPL [keywd=value ...]

where stub was specified in

        ampl -obstub ...
or
        ampl -ogstub...

Such an invocation causes baron to read from stub.nl and to write stub.sol.

-----------------
Controlling baron
-----------------

Baron reads keywords and values from the environment (shell) variable
baron_options and from the command line.  Execute

        baron -?

or (if your shell requires ? to be quoted)

        baron '-?'

for a summary of baron usage and

        baron -=

(or baron '-=') for a summary of keywords peculiar to the AMPL/BARON
driver "baron".

------------------
Sample Invocations
------------------

  If you're using AMPL, just say

        option solver baron;
        solve;

  If you've executed, say,

        ampl -objunk junk.model junk.data

then you could say

        baron junk epsr=1e-4 lsolver=minos

to force baron to use relative tolerance .0001 between the function
value at the solution returned and a lower bound thereon and to use
local solver minos.  (This assumes a license for the AMPL/MINOS solver
"minos" is available.)  With the Bourne shell, either of the invocations

        baron_options='epsr=1e-4 lsolver=minos' baron junk
or
        baron_options='epsr=1e-4 lsolver=minos'
        export baron_options
        baron junk

would have the same effect; within AMPL, specifying

        option baron_options 'epsr=1e-4 lsolver=minos', solver baron;
        solve;

would also have this effect.


-----------------------
solve_result_num values
=======================

Here is a table of solve_result_num values that "baron" can return to
an AMPL session, along with the text that appears in the associated
solve_message.

        Value   Message

        0       optimal within tolerances
        100     numerical difficulties but possibly optimal
	150	feasible but insufficient progress
        200     infeasible
        201     infeasible, IIS found
        202     infeasible, IS found, possibly not irreducible
        203     infeasible, IIS sought but not found
        300     unbounded
        400     node limit reached
        401     iteration limit reached
        402     CPU time limit reached
        500     licensing error
        501     numerical difficulties
        502     interrupted (Control-C)
        503     too little memory
        504     terminated by BARON
        505     BARON syntax error (should not happen)
        506     operation not supported by BARON
        507     Interrupted by Control-C

----------
INSTALLING
==========

On Linux systems, libbaron*.so (where the values of "*" depends on the
current version of BARON) needs to appear in the current directory
when baron itself appears there, or in one of the standard places
(specified by /etc/ld.so.conf on some systems), or in a directory
named in $LD_LIBRARY_PATH.  An alternative is to add a short shell
script, such as

        #!/bin/sh
        LD_LIBRARY_PATH=/usr/local/lib
        export LD_LIBRARY_PATH
        exec /usr/local/bin/baronx "$@"

to a directory in your usual $PATH (and mark the script executable
with, e.g., "chmod +x baron").  The above script assumes that the
true "baron" binary has been moved to /usr/local/bin/baronx and that
the libbaron*.so file has been moved to /usr/local/lib.

MacOSX systems are similar to Linux systems, but with libbaron*.dylib
in place of libbaron*.so and with DYLD_LIBRARY_PATH in place of
LD_LIBRARY_PATH.  On MacOSX systems, it suffices for the
libbaron*.dylib file to appear in the same directory as baron, at
least when that directory is current when baron is invoked.

On MS Windows systems, xpress.exe and the relevant baron-*.dll must
appear somewhere in your usual search $PATH (or in the current
directory).

Questions about this stuff?  Contact support@ampl.com .

