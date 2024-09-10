#include "mp2nlmodelapi.h"
#include "mp/nl-solver.hpp"
#include "mp/nl-opcodes.h"


namespace mp {

void MP2NLModelAPI::InitProblemModificationPhase(const FlatModelInfo* flat_model_info) {
  alg_con_info_.reserve(
      flat_model_info->GetNumberOfConstraintsOfGroup(CG_Algebraic));
  log_con_info_.reserve(
      flat_model_info->GetNumberOfConstraintsOfGroup(CG_Logical));
  sos_info_.reserve(
      flat_model_info->GetNumberOfConstraintsOfGroup(CG_SOS));
}

void MP2NLModelAPI::AddVariables(const VarArrayDef& vad) {
  var_lbs_ = {vad.plb(), (size_t)vad.size()};
  var_ubs_ = {vad.pub(), (size_t)vad.size()};
  var_types_ = {vad.ptype(), (size_t)vad.size()};
  var_names_ = {vad.pnames(), (size_t)vad.size()};
  mark_data_.col_sizes_orig_.resize(var_lbs_.size());
}


void MP2NLModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  assert(iobj == (int)obj_info_.size());
  obj_info_.push_back(MakeItemInfo(lo, StaticItemTypeID::ID_LinearObjective));
}

void MP2NLModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  assert(iobj == (int)obj_info_.size());
  /// @todo ?
  throw std::runtime_error("Quadratic objective not supported");
}

void MP2NLModelAPI::SetNLObjective( int iobj, const NLObjective& nlo ) {
  assert(iobj == (int)obj_info_.size());
  obj_info_.push_back(MakeItemInfo(nlo, StaticItemTypeID::ID_NLObjective));
}


MP2NL_Expr MP2NLModelAPI::GetVarExpression(int i) {
  return {i+1};    // ?
}

MP2NL_Expr MP2NLModelAPI::GetZeroExpression() {
  return {};
}

void MP2NLModelAPI::AddConstraint(const LinConRange& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConRange)); }

void MP2NLModelAPI::AddConstraint(const LinConLE& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConLE)); }

void MP2NLModelAPI::AddConstraint(const LinConEQ& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConEQ)); }

void MP2NLModelAPI::AddConstraint(const LinConGE& lc)
{ alg_con_info_.push_back(MakeItemInfo(lc, StaticItemTypeID::ID_LinConGE)); }


/// To access information from an NLConstraint,
/// use the following accessors (don't use methods of NLConstraint itself):
/// - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
///   GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
///
/// Implementation follows partly reader_nl.cc from SCIP.
void MP2NLModelAPI::AddConstraint( const NLConstraint& nlc )
{ alg_con_info_.push_back(MakeItemInfo(nlc, StaticItemTypeID::ID_NLConstraint)); }

void MP2NLModelAPI::AddConstraint( const NLAssignEQ& nlae )
{ alg_con_info_.push_back(MakeItemInfo(nlae, StaticItemTypeID::ID_NLAssignEQ)); }
void MP2NLModelAPI::AddConstraint( const NLAssignLE& nlae )
{ alg_con_info_.push_back(MakeItemInfo(nlae, StaticItemTypeID::ID_NLAssignLE)); }
void MP2NLModelAPI::AddConstraint( const NLAssignGE& nlae )
{ alg_con_info_.push_back(MakeItemInfo(nlae, StaticItemTypeID::ID_NLAssignGE)); }

void MP2NLModelAPI::AddConstraint(const NLComplementarity& cc)
{ alg_con_info_.push_back(MakeItemInfo(cc, StaticItemTypeID::ID_NLComplementarity)); }


void MP2NLModelAPI::AddConstraint( const NLLogical& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLLogical)); }

void MP2NLModelAPI::AddConstraint( const NLEquivalence& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLEquivalence)); }
void MP2NLModelAPI::AddConstraint( const NLImpl& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLImpl)); }
void MP2NLModelAPI::AddConstraint( const NLRimpl& nll )
{ log_con_info_.push_back(MakeItemInfo(nll, StaticItemTypeID::ID_NLRimpl)); }


void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinLE &ic)
{ log_con_info_.push_back(MakeItemInfo(ic, StaticItemTypeID::ID_IndicatorConstraintLinLE)); }
void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinEQ &ic)
{ log_con_info_.push_back(MakeItemInfo(ic, StaticItemTypeID::ID_IndicatorConstraintLinEQ)); }
void MP2NLModelAPI::AddConstraint(const IndicatorConstraintLinGE &ic)
{ log_con_info_.push_back(MakeItemInfo(ic, StaticItemTypeID::ID_IndicatorConstraintLinGE)); }



void MP2NLModelAPI::AddConstraint(const SOS1Constraint& sos)
{ sos_info_.push_back(MakeItemInfo(sos, StaticItemTypeID::ID_SOS1Constraint)); }
void MP2NLModelAPI::AddConstraint(const SOS2Constraint& sos)
{ sos_info_.push_back(MakeItemInfo(sos, StaticItemTypeID::ID_SOS2Constraint)); }


MP2NL_Expr MP2NLModelAPI::AddExpression(const LinExpression &le) {
  return {};
}

MP2NL_Expr MP2NLModelAPI::AddExpression(const QuadExpression &qe) {
  return {};
}


MP2NL_Expr MP2NLModelAPI::AddExpression(const AbsExpression &abse) {
  return {};
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




NLHeader MP2NLModelAPI::Header() {
  PrepareModel();
  return DoMakeHeader();
}

void MP2NLModelAPI::PrepareModel() {
  MarkVars();
  SortVars();
  MarkItems();
}

void MP2NLModelAPI::MarkVars() {
  mark_data_.var_prior_.resize(var_lbs_.size());
  for (auto i=var_lbs_.size(); i--; ) {
    mark_data_.var_prior_[i].second = i;
    if (var::Type::INTEGER == var_types_[i]) {
      if (!var_lbs_[i] && 1.0==var_ubs_[i]) {
        mark_data_.var_prior_[i].first = 1;
        ++mark_data_.n_var_lin_bin_;
      } else {
        mark_data_.var_prior_[i].first = 2;
        ++mark_data_.n_var_lin_int_;
      }
    }
  }
}

void MP2NLModelAPI::SortVars() {
  std::sort(mark_data_.var_prior_.begin(), mark_data_.var_prior_.end());
  mark_data_.var_order_12_.resize(var_lbs_.size());
  mark_data_.var_order_21_.resize(var_lbs_.size());
  for (auto i=var_lbs_.size(); i--; ) {
    mark_data_.var_order_12_[i] = mark_data_.var_prior_[i].second;
    mark_data_.var_order_21_[mark_data_.var_prior_[i].second] = i;
  }
}

void MP2NLModelAPI::MarkItems() {
  ItemMarkingData prm;
  for (auto& info: alg_con_info_)
    info.GetDispatcher().MarkItem(info.GetPItem(), prm);
}

void MP2NLModelAPI::Add2ColSizes(const std::vector<int>& vars) {
  for (auto v: vars)
    ++mark_data_.col_sizes_orig_[v];
}


NLHeader MP2NLModelAPI::DoMakeHeader() {
  mp::NLHeader hdr;

  hdr.num_vars = var_lbs_.size();
  hdr.num_algebraic_cons = alg_con_info_.size();
  hdr.num_objs = obj_info_.size();
  hdr.num_ranges = mark_data_.n_ranges_;
  hdr.num_eqns = mark_data_.n_eqns_;
  hdr.num_logical_cons = 0;

  /** Total number of nonlinear constraints. */
  hdr.num_nl_cons = 0;
  hdr.num_nl_objs = 0;
  hdr.num_compl_conds = 0;
  hdr.num_nl_compl_conds = 0;
  hdr.num_compl_dbl_ineqs = 0;
  hdr.num_compl_vars_with_nz_lb = 0;

  /** Number of nonlinear network constraints. */
  hdr.num_nl_net_cons = 0;
  hdr.num_linear_net_cons = 0;

  /**
      Number of nonlinear variables in constraints including nonlinear
      variables in both constraints and objectives.
     */
  hdr.num_nl_vars_in_cons = 0;

  /**
      Number of nonlinear variables in objectives including nonlinear
      variables in both constraints and objectives.
     */
  hdr.num_nl_vars_in_objs = 0;

  /** Number of nonlinear variables in both constraints and objectives. */
  hdr.num_nl_vars_in_both = 0;

  // Miscellaneous
  // -------------

  /** Number of linear network variables (arcs). */
  hdr.num_linear_net_vars = 0;

  /** Number of functions. */
  hdr.num_funcs = 0;

  // Information about discrete variables
  // ------------------------------------

  /** Number of linear binary variables. */
  hdr.num_linear_binary_vars = mark_data_.n_var_lin_bin_;

  /** Number of linear non-binary integer variables. */
  hdr.num_linear_integer_vars = mark_data_.n_var_lin_int_;

  /**
      Number of integer nonlinear variables in both constraints and objectives.
     */
  hdr.num_nl_integer_vars_in_both = 0;

  /** Number of integer nonlinear variables just in constraints. */
  hdr.num_nl_integer_vars_in_cons = 0;

  /** Number of integer nonlinear variables just in objectives. */
  hdr.num_nl_integer_vars_in_objs = 0;

  // Information about nonzeros
  // --------------------------

  /** Number of nonzeros in constraints' Jacobian. */
  hdr.num_con_nonzeros = 0;

  /** Number of nonzeros in all objective gradients. */
  hdr.num_obj_nonzeros = 0;

  // Information about names
  // -----------------------

  /** Length of longest con/obj name if names are available. */
  hdr.max_con_name_len = 0;    // no need to set

  /** Length of longest variable name if names are available. */
  hdr.max_var_name_len = 0;    // no need to set

  // Information about common expressions
  // ------------------------------------

  /**
      Number of common expressions that appear both in constraints
      and objectives.
     */
  hdr.num_common_exprs_in_both = 0;

  /**
      Number of common expressions that appear in multiple constraints
      and don't appear in objectives.
     */
  hdr.num_common_exprs_in_cons = 0;

  /**
      Number of common expressions that appear in multiple objectives
      and don't appear in constraints.
     */
  hdr.num_common_exprs_in_objs = 0;

  /**
      Number of common expressions that only appear in a single constraint
      and don't appear in objectives.
     */
  hdr.num_common_exprs_in_single_cons = 0;

  /**
      Number of common expressions that only appear in a single objective
      and don't appear in constraints.
     */
  hdr.num_common_exprs_in_single_objs = 0;

  hdr.prob_name = "mp2nl_model";

  return hdr;
}


void CreateNLS() {
  Env env;
  MP2NLModelAPI nlf(env);
  SOLHandler esolh;
  mp::NLUtils utils;

  mp::NLSolver nlsol(&utils);

  std::string solver = "minos";
  std::string sopts = "";

  if (!nlsol.Solve(nlf, esolh, solver, sopts)) {
    printf("%s\n", nlsol.GetErrorMessage());
  } else {
  }
}



} // namespace mp
