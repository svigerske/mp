################################################
## defvar_01.mod: test def var elimination,
## depending on usage.
## See sp/#81.
################################################

var x <=0;
var y >=0;
var dx = 0.1*x + 2.4;
var dy = 0.2*y + 8;
var dxy = 0.8 + (2.6 - dx)*(dy*1.2);

s.t. L1: y - x >= 3;

minimize O1: dxy;
