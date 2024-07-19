##########################################
## Test products of binary variables, in particular reformulations into logicals
## prod_bin_03.mod
##########################################

var b {1..3} binary;

minimize Obj: b[1] - b[2] + b[3];

s.t. C1: b[1]*(1-b[2])*b[3] == 1;
