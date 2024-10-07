####################################
# Test that we don't convert abs(x) <= 5 into a cone
# TODO Also perform the actual test, now it's just for solving
####################################

var b1 binary;
var b2 binary;
var x >=-19 <=32;
var y >= -12 <= 5;

minimize Obj: b1 - b2 - x + 3*y;

var charge_diff = abs(x-y);

s.t. Con: b1 && b2 ==> x<=y;

s.t. Cone: charge_diff <= if b1>=1 then 5 else 15;
