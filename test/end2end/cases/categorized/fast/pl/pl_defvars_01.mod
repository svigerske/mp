## pl_defvars_01.mod
## Using PL terms in defined variables

var x {1..3} >= -5 <= 23.4;

var f1 = <<3,8; 1,0.1,-2>> x[1];

var f3 = <<-1,4.7,12; -6,5,-4,7.9>> x[3];

var f2 = <<-2, 6.9; -1.3, -.03, 1.5>> x[2];

minimize O1: 5*f1 + 4*f2;

s.t. LC1: f1+f2 - f3 >= 0.78;
