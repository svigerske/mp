#ifndef MP2NLCOMMON_H
#define MP2NLCOMMON_H

#include <vector>
#include <string>

#include "mp/backend-to-model-api.h"

#include "mp/format.h"


namespace mp {

/// Information shared by both
/// `ScipBackend` and `ScipModelAPI`
struct MP2NLCommonInfo {

private:
};


/// Common API for MP2NL classes
class MP2NLCommon :
		public Backend2ModelAPIConnector<MP2NLCommonInfo> {
public:
  static constexpr double Infinity() { return INFINITY; }
  static constexpr double MinusInfinity() { return -INFINITY; }

protected:
  void OpenSolver();
  void CloseSolver();

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;
};

} // namespace mp

#endif // MP2NLCOMMON_H
