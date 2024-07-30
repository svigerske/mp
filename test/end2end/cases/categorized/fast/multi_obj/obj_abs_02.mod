#######################################
## obj_abs_02.mod
## Test objective sense-aware reformulation
## Test .objweight in obj:multi mode
## abs()
#######################################

var x >=-2 <=5;
var y >=-12 <=3;


suffix objpriority IN;
suffix objweight IN;

## The 1st objective is dummy, to enable the multi-objective approach
minimize ObjD: 3*x - 2*y suffix objpriority 1, suffix objweight -1;
maximize Obj1: abs(x) - 2*abs(y) suffix objpriority 2, suffix objweight 3;

s.t. C1: 3*x - 2*y == 13;

