#include "mp/format.h"
#include "mp2nlcommon.h"


namespace mp {

void MP2NLCommon::OpenSolver() { }

void MP2NLCommon::CloseSolver() { }

int MP2NLCommon::NumLinCons() const {
  return 0;
}

int MP2NLCommon::NumVars() const {
  return 0.0;
}

int MP2NLCommon::NumObjs() const {
  return 1;
}

int MP2NLCommon::NumQPCons() const {
  int count = 0;
  return count;
}

int MP2NLCommon::NumSOSCons() const {
  int count = 0;
  return count;
}

int MP2NLCommon::NumIndicatorCons() const {
  int count = 0;
  return count;
}

} // namespace mp
