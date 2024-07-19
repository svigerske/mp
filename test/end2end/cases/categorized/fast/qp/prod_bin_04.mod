##########################################
## Test products of binary variables, in particular reformulations into logicals
## prod_bin_04.mod
##########################################

var b {1..5} binary;

minimize Obj: b[1] - b[2] + b[3] - b[4] + b[5];

s.t. C1: b[1]*((b[2]-1)*(b[3]-1)-1)*b[4] == -1;

s.t. C2: b[1]*(1-b[5])*(b[4]+15) >= 7;