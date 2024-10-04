
# -------------------------------------------------------------
# x!=const for float variable with 2 constants
# -------------------------------------------------------------

var x >= -2e3, <= 1e3;
var y;

minimize XMinus5:
    y;

subj to NE1:
    x!=(-1e2-2);

subj to NE2:
    x!=-100;

subj to YUp:
    y >= x-(-1e2-2);

subj to YDown:
    y >= (-1e2-2)-x;
