#### See if Div/const is well presolved
#### MP does this because AMPL 20240331 does not
#### See if abs() is well redefined (or used well natively)
###################################################

var x {1..3} >=-200, <=200;

minimize Unbalance:
   sum {i in 1..3} abs(x[i] - sum {i1 in 1..3} x[i1] / 3);

## Some extra constraints
 s.t. C1: sum {I in 1..3} x[I] == 15;

 s.t. C2: x[1] + 2*x[2] + 4*x[3] == 18;

