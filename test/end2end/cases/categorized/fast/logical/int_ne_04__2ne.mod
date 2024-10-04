##############################
## Conditional disequality
##############################

var x {1..2} >=-5 <=18 integer;
var y {1..2} binary;

minimize TotalDiff: abs(x[1]-x[2]) + abs(y[1]-y[2]);

s.t. OneNE: x[1]!=x[2] || y[1]!=y[2];
