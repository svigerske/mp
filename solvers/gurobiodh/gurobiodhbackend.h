#ifndef MP_GUROBIODH_BACKEND_H_
#define MP_GUROBIODH_BACKEND_H_

#include <string>
#include "odh/odhcommon.h"
#include "gurobi/gurobibackend.h"


namespace mp {

  class GurobiODHBackend : public GurobiBackend, ODHCommonInfo {
  
    /// Chance for the Backend to init solver environment, etc
    void InitOptionParsing() override;

    void Solve() override;

    ~GurobiODHBackend();
  
  };
}
#endif