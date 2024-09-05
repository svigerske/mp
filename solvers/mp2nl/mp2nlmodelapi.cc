#include "mp2nlmodelapi.h"


namespace mp {

void MP2NLModelAPI::InitProblemModificationPhase(const FlatModelInfo* flat_model_info) {
}

void MP2NLModelAPI::AddVariables(const VarArrayDef& v) {
}

void MP2NLModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
}


void MP2NLModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  /// @todo
  throw std::runtime_error("Quadratic objective not supported");
}


MP2NL_Expr MP2NLModelAPI::GetVarExpression(int i) {
  return {i+1};    // ?
}

MP2NL_Expr MP2NLModelAPI::GetZeroExpression() {
  return {};
}

void MP2NLModelAPI::AddConstraint(const LinConRange& lc) {
}

void MP2NLModelAPI::AddConstraint(const LinConLE& lc) {
}

void MP2NLModelAPI::AddConstraint(const LinConEQ& lc) {
}

void MP2NLModelAPI::AddConstraint(const LinConGE& lc) {
}


MP2NL_Expr MP2NLModelAPI::AddExpression(const AbsExpression &abse) {
  return {};
}


void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)  {
}

void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)  {
}

void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)  {
}


/// To access information from an NLConstraint,
/// use the following accessors (don't use methods of NLConstraint itself):
/// - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
///   GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
///
/// Implementation follows partly reader_nl.cc from SCIP.
void MP2NLModelAPI::AddConstraint( const NLConstraint& nlc ) {
}

void MP2NLModelAPI::AddConstraint( const NLAssignEQ& nlae ) {
}
void MP2NLModelAPI::AddConstraint( const NLAssignLE& nlae ) {
}
void MP2NLModelAPI::AddConstraint( const NLAssignGE& nlae ) {
}

void MP2NLModelAPI::AddConstraint( const NLLogical& nll ) {
}

void MP2NLModelAPI::AddConstraint( const NLEquivalence& nll ) {
}
void MP2NLModelAPI::AddConstraint( const NLImpl& nll ) {
}
void MP2NLModelAPI::AddConstraint( const NLRimpl& nll ) {
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const LinExpression &le) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const QuadExpression &qe) {
  return {};
}


void MP2NLModelAPI::AddConstraint(const SOS1Constraint& sos) {
}

void MP2NLModelAPI::AddConstraint(const SOS2Constraint& sos) {
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const AndExpression &ee) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const OrExpression &ee) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const ExpExpression &ee) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const LogExpression &ee) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const PowExpression &ee) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const SinExpression &ee) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const CosExpression &ee) {
  return {};
}

void MP2NLModelAPI::FinishProblemModificationPhase() {
}


} // namespace mp
