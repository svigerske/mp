# expcones_01__plain.mod
# From https://docs.mosek.com/latest/capi/tutorial-ceo-shared.html

var x {i in 1..3} >= if i<3 then 0 else -Infinity;

minimize Obj:
   x[1] + x[2];

s.t. ExpCone:
   x[1] >= x[2] * exp( x[3] / x[2] );

s.t. LinCon:
   sum {i in 1..3} x[i] == 1;
