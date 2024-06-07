# expcones_03.mod

var x {i in 1..3} >= if i<3 then 0 else -Infinity;

minimize Obj:
   x[1] + x[2];

s.t. ExpCone:
   2*x[1] >= 3*x[2] * exp( 4*x[3] / (3*x[2]) );

s.t. LinCon:
   sum {i in 1..3} x[i] == 1;
