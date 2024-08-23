#include "cplex/cplex-ampls-c-api.h"

AMPLS_C_EXPORT AMPLS_MP_Solver* AMPLSOpen_cplex(int argc, char** argv)
{
  CCallbacks cb = { NULL };
  return Open_cplex(cb);
}