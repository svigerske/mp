#include "mp/format.h"
#include "cplexcommon.h"

namespace mp {


int CplexCommon::NumLinCons() const {
  return CPXgetnumrows (env(), lp());
}

int CplexCommon::NumVars() const {
  return CPXgetnumcols (env(), lp());
}

int CplexCommon::NumObjs() const {
  return CPXgetnumobjs (env(), lp());
}

int CplexCommon::NumQPCons() const {
  return CPXgetnumqconstrs(env(), lp());
}

int CplexCommon::NumIndicatorCons() const {
  return CPXgetnumindconstrs(env(), lp());
}

int CplexCommon::NumSOSCons() const {
  return CPXgetnumsos(env(), lp());
}
bool CplexCommon::HasQObj() const {
  return CPXgetnumquad(env(), lp()) > 0;
}
int CplexCommon::ModelSense() const {
  return CPXgetobjsen(env(), lp());
}

void CplexCommon::GetSolverOption(int key, int &value) const {
  CPLEX_CALL( CPXgetintparam(env(), key, &value) );
}

void CplexCommon::SetSolverOption(int key, int value) {
  CPLEX_CALL( CPXsetintparam(env(), key, value) );
}

void CplexCommon::GetSolverOption(int key, double &value) const {
  CPLEX_CALL( CPXgetdblparam(env(), key, &value) );
}

void CplexCommon::SetSolverOption(int key, double value) {
  CPLEX_CALL( CPXsetdblparam(env(), key, value) );
}

void CplexCommon::GetSolverOption(int key, std::string &value) const {
  char buffer[CPX_STR_PARAM_MAX];
  CPLEX_CALL( CPXgetstrparam(env(), key, buffer) );
  value = buffer;
}

void CplexCommon::SetSolverOption(int key, const std::string& value) {
  CPLEX_CALL( CPXsetstrparam(env(), key, value.c_str()) );
}

double CplexCommon::GetCPLEXDblParam(int param) {
  double value;
  CPLEX_CALL(CPXgetdblparam(env(), param, &value));
  return value;
}
int CplexCommon::GetCPLEXIntParam(int param) {
  int value;
  CPLEX_CALL(CPXgetintparam(env(), param, &value));
  return value;
}
void CplexCommon::SetCPLEXParam(int param, int value) {
  CPLEX_CALL(CPXsetintparam(env(), param, value));
}
void CplexCommon::SetCPLEXParam(int param, double value) {
  CPLEX_CALL(CPXsetdblparam(env(), param, value));
}
} // namespace mp
