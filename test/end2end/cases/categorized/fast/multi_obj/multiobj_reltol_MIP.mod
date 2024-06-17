
var x;
var y;
var z integer;       # For Gurobi to handle AbsTol in direct way

maximize xobj: x;
minimize yobj: y;

subj to c1: y >= x;
subj to c2: y <= x + 1;

subj to c3: y <= 2*x;
subj to c4: y >= 2*x - 2;

subj to c5_keep_z: y >= z;

suffix objpriority IN;
suffix objreltol IN;

let xobj.objpriority := 10;
let yobj.objpriority := 1;

let xobj.objreltol:= 0.25/3;

# Need mip:gap=0 due to the way tolerances are applied:
# https://www.gurobi.com/documentation/current/refman/working_with_multiple_obje.html.
# option gurobi_options "obj:multi=1 writeprob=multiobj1_MIP.lp mip:gap=0";
# option solver gurobi;

# solve;
# display x, y, z;
