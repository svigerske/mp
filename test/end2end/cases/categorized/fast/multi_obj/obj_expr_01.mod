## Nonlinear expressions in objectives

var x {1..3} >=-30 <=54;

s.t. X1LB: x[1] >=34;
s.t. X2LB: x[2] >= 2;

s.t. ConL1L: x[1] + 2*x[2] <= 137;

s.t. ConNL1: x[1] + x[3] - x[2]^1.5 >= -17;

suffix objpriority;          ## To enable lexicographic optimization

maximize Obj1: -0.55*x[1] + 4*x[2] + log(abs(x[3] / x[2]) + 0.5)
   suffix objpriority 2;

minimize Obj4: 5*x[2] + x[3] + x[2]^1.24
   suffix objpriority 1;
