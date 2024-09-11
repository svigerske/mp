#include "mp2nlmodelapi.h"
#include "mp/nl-solver.hpp"
#include "mp/nl-opcodes.h"
#include "mp/sol-handler.h"


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


/// Implement NLSolverIntf
class MP2NLSolverImpl
    : public MP2NLSolverIntf,
      public SOLHandler {
public:
  /// Construct
  MP2NLSolverImpl(MP2NLModelAPI& mapi) : mapi_(mapi) { }

  /// Solve
  void Solve(const char* solver, const char* sopts) override {
    if_solve_attempted_ = true;
    if_solve_ok_ = nlsol_.Solve(mapi_, *this, solver, sopts);
  }

  /// AMPL solve result code
  int GetSolveResult() const override {
    return (if_solve_ok_ || !if_solve_attempted_)
        ? sresult_ : 500;             // 500: failure
  }

  /// Solve result message or error message
  const char* GetSolveMessage() const override {
    if (if_solve_ok_)
      solve_message_final_ = solve_message_;
    else {
      solve_message_final_ = nlsol_.GetErrorMessage();
      if (solve_message_.size()) {
        solve_message_final_ +=
            "\nOriginal solve message: ";
        solve_message_final_ += solve_message_;
      }
    }
    return solve_message_final_.c_str();
  }

  /// Number of backspaces to print
  /// if printing the solve message right here,
  /// or skip so many symbols first.
  int GetSolveMessageNbs() const override { return nbs_; }

  /// Stub file used
  const char* GetFileStub() const override
  { return nlsol_.GetFileStub().c_str(); }

  /// Objno used
  int GetObjnoUsed() const override { return objno_used_; }

  /// Primal solution
  ArrayRef<double> GetX() const override { return primals_; }

  /// Dual solution
  ArrayRef<double> GetY() const override { return duals_; }

  /// @todo + Suffixes... Pull or push?
  /// ..

  /////////////////////////////////////////////////////////////////////////
  //////////////////////// SOLHandler implementation //////////////////////
  /////////////////////////////////////////////////////////////////////////
public:
  /** The NLHeader used to write the NL file. */
  NLHeader Header() const { return mapi_.Header(); }

  /** Receive solve message.
   *  The message always ends with '\n'.
   *
   *  @param nbs: number of backspaces
   *  in the original solve message.
   *  So many characters should be skipped
   *  from the message if printed straightaway.
   *  AMPL solver drivers can supply the message
   *  with initial backspaces to indicate
   *  that so many characters should be skipped
   *  when printing. For example, if the driver prints
   *  MINOS 5.51:
   *  and exits, and the message starts with that again,
   *  this part should be skipped.
   */
  void OnSolveMessage(const char* s, int nbs) {
    solve_message_ = s;
    nbs_ = nbs;
  }

  /**
   * Can be ignored by external systems.
   * @return non-zero to stop solution input.
   */
  int OnAMPLOptions(const AMPLOptions& ) { return 0; }

  /**
   * Dual values for algebraic constraints,
   * if provided in the solution.
   * Number of values <= NumAlgCons().
   * Implementation:
   *
   *   duals.reserve(rd.Size());
   *   while (rd.Size())
   *     duals.push_back(rd.ReadNext());
   */
  template <class VecReader>
  void OnDualSolution(VecReader& rd) {
    duals_.clear();
    duals_.reserve(rd.Size());
    while (rd.Size())
      duals_.push_back( rd.ReadNext() );
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& rd) {
    primals_.clear();
    primals_.reserve(rd.Size());
    while (rd.Size())
      primals_.push_back( rd.ReadNext() );
  }

  /**
   * Receive notification of the objective index
   * used by the driver (solver option 'objno'-1).
   */
  void OnObjno(int objno) { objno_used_ = objno; }

  /**
   * Receive notification of the solve code.
   * Solve result codes docu:
   * https://mp.ampl.com/features-guide.html#solve-result-codes
   */
  void OnSolveCode(int sr) { sresult_ = sr; }

  /**
   * OnIntSuffix().
   *
   * For constraints, can include values for
   * logical constraints (after algebraic.)
   * Sparse representation - can be empty
   * (i.e., all values zero.)
   *
   * const auto& si = sr.SufInfo();
   * int kind = si.Kind();
   * int nmax = nitems_max[kind & 3];
   * const std::string& name = si.Name();
   * const std::string& table = si.Table();
   * while (sr.Size()) {
   *   std::pair<int, int> val = sr.ReadNext();
   *   if (val.first<0 || val.first>=nmax) {
   *     sr.SetError(NLW2_SOLRead_Bad_Suffix,
   *       "bad suffix element index");
   *     return;
   *   }
   *   suf[val.first] = val.second;
   * }
   * if (NLW2_SOLRead_OK == sr.ReadResult())    // Can check
   *   RegisterSuffix(kind, name, table, suf);
   */
  template <class SuffixReader>
  void OnIntSuffix(SuffixReader& sr) {
    while (sr.Size())
      sr.ReadNext();       // read & forget by default
  }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
  template <class SuffixReader>
  void OnDblSuffix(SuffixReader& sr) {
    while (sr.Size())
      sr.ReadNext();       // read & forget by default
  }


private:
  MP2NLModelAPI& mapi_;

  mp::NLUtils utils_;

  mp::NLSolver nlsol_ { &utils_ };


  /// Solution
  bool if_solve_attempted_ {};
  bool if_solve_ok_ {};     // return NLSolver's error message if not
  std::string solve_message_;
  mutable std::string solve_message_final_;
  int nbs_ {};
  std::vector<double> duals_,
      primals_;
  int objno_used_ {-1},
      sresult_ {-1};
};


void MP2NLModelAPI::CreateInterfaces() {
  p_nls_ = std::make_unique<MP2NLSolverImpl>(*this);
  this->MP2NLCommonInfo::p_nlsi_ = p_nls_.get();
}


} // namespace mp
