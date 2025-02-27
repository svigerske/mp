#include "coptmodelapi.h"

namespace mp {

void CoptModelAPI::InitProblemModificationPhase(const FlatModelInfo*) { }

void CoptModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<char> vtypes(v.size());
  for (size_t i=v.size(); i--; )
    vtypes[i] = var::Type::CONTINUOUS==v.ptype()[i] ?
          COPT_CONTINUOUS : COPT_INTEGER;
  COPT_CCALL(COPT_AddCols(lp(), (int)v.size(), NULL, NULL,
    NULL, NULL, NULL, vtypes.data(), v.plb(), v.pub(),  v.pnames()));
}

void CoptModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj<1) {
    COPT_CCALL(COPT_SetObjSense(lp(), 
                    obj::Type::MAX==lo.obj_sense() ? COPT_MAXIMIZE : COPT_MINIMIZE) );
    double zero_out = 0.0;
    for (int i=NumVars(); i--; )
      COPT_CCALL(COPT_SetColObj(lp(), 1, &i, &zero_out));
    COPT_CCALL(COPT_SetColObj(lp(), lo.num_terms(),
                           lo.vars().data(), lo.coefs().data()) );
  } else {
//    TODO
  }
}


void CoptModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);                         // add the linear part
    const auto& qt = qo.GetQPTerms();
    COPT_CCALL(COPT_SetQuadObj(lp(), qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
      (double*)qt.pcoefs()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void CoptModelAPI::AddConstraint(const LinConRange& lc) {
	COPT_CCALL(COPT_AddRow(lp(), lc.size(),
												 lc.pvars(), lc.pcoefs(), 0,
												 lc.lb() < -COPT_INFINITY ? -COPT_INFINITY : lc.lb(),
												 lc.ub() > COPT_INFINITY ? COPT_INFINITY : lc.ub(),
												 lc.name()));
}

void CoptModelAPI::AddConstraint(const LinConLE& lc) {
  char sense = COPT_LESS_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
		sense, lc.rhs(), 0, lc.name()));
}
void CoptModelAPI::AddConstraint(const LinConEQ& lc) {
  char sense = COPT_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
		sense, lc.rhs(), 0, lc.name()));
}
void CoptModelAPI::AddConstraint(const LinConGE& lc) {
  char sense = COPT_GREATER_EQUAL;
  COPT_CCALL(COPT_AddRow(lp(), lc.size(), lc.pvars(), lc.pcoefs(),
		sense, lc.rhs(), 0, lc.name()));
}

void CoptModelAPI::AddConstraint(const QuadConLE& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  COPT_CCALL(COPT_AddQConstr(lp(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                             qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
														 (double*)qt.pcoefs(), COPT_LESS_EQUAL, qc.rhs(), qc.name()));
}

void CoptModelAPI::AddConstraint(const QuadConEQ& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  COPT_CCALL(COPT_AddQConstr(lp(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                             qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
														 (double*)qt.pcoefs(), COPT_EQUAL, qc.rhs(), qc.name()));
}

void CoptModelAPI::AddConstraint(const QuadConGE& qc) {
  const auto& lt = qc.GetLinTerms();
  const auto& qt = qc.GetQPTerms();
  COPT_CCALL(COPT_AddQConstr(lp(), lt.size(), (int*)lt.pvars(), (double*)lt.pcoefs(),
                             qt.size(), (int*)qt.pvars1(), (int*)qt.pvars2(),
														 (double*)qt.pcoefs(), COPT_GREATER_EQUAL, qc.rhs(), qc.name()));
}



void CoptModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
  COPT_CCALL(COPT_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    COPT_LESS_EQUAL,
    ic.get_constraint().rhs()));
                               
}
void CoptModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
  COPT_CCALL(COPT_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    COPT_EQUAL,
    ic.get_constraint().rhs()));
}
void CoptModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
  COPT_CCALL(COPT_AddIndicator(lp(),
    ic.get_binary_var(), ic.get_binary_value(),
    (int)ic.get_constraint().size(),
    ic.get_constraint().pvars(),
    ic.get_constraint().pcoefs(),
    COPT_GREATER_EQUAL,
    ic.get_constraint().rhs()));
}

void CoptModelAPI::AddConstraint(const SOS1Constraint& sos) {
  int type = COPT_SOS_TYPE1;
  int beg = 0;
  const int size = sos.size();
  COPT_CCALL(COPT_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));
}

void CoptModelAPI::AddConstraint(const SOS2Constraint& sos) {
  int type = COPT_SOS_TYPE2;
  int beg = 0;
  const int size = sos.size();
  COPT_CCALL(COPT_AddSOSs(lp(), 1, &type, &beg,
    &size, (int*)sos.get_vars().data(),
    (double*)sos.get_weights().data()));
}


void CoptModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
