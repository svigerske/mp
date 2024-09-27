#include "mp2nlmodelapi.h"
#include "mp/nl-solver.hpp"
#include "mp/nl-opcodes.h"
#include "mp/sol-handler.h"


namespace mp {

void MP2NLModelAPI::InitCustomOptions() {
  GetEnv().AddStoredOption("nl:stub stub nlstub",
                           "Filename stub for the NL and SOL files.",
                           storedOptions_.stub_);
  GetEnv().AddStoredOption("nl:text nltextformat nltext",
                           "0*/1: write NL file in text format.",
                           storedOptions_.nl_format_text_);
  GetEnv().AddStoredOption("nl:comments nlcomments",
                           "0*/1: add comments to the text-format NL file.",
                           storedOptions_.nl_comments_);

}

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
  if (vad.pnames())
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


template <class Expr>
MP2NL_Expr MP2NLModelAPI::AddExpression(
    const Expr &expr, ExpressionTypeID eid) {
  // Only should visit such arguments
  // which are true expressions, i.e.,
  // not the pure-variable parts.
  // Thus, need to specialize for NLDefVar.
  VisitArguments(expr, [this](MP2NL_Expr mp2nle){
    RegisterExpression(mp2nle);
  });
  return StoreMP2NLExprID(expr, eid);
}

template <class Expr>
MP2NL_Expr MP2NLModelAPI::StoreMP2NLExprID(
    const Expr &expr, ExpressionTypeID eid) {
  expr_info_.push_back(MakeItemInfo(expr, eid));
  expr_counter_.push_back(0);
  return MakeExprID( int(expr_info_.size()-1) );
}

void MP2NLModelAPI::RegisterExpression(MP2NL_Expr expr) {
  CountExpression(expr);
}

/// Count expression depending on its kind.
/// This duplicates the counters in FlatConverter
///   -- could check equality.
void MP2NLModelAPI::CountExpression(MP2NL_Expr expr) {
  if (expr.IsExpression()) {
    auto index = expr.GetExprIndex();
    assert(index >=0 && index < expr_counter_.size()
           && index < expr_info_.size());
    ++expr_counter_[index];
  }   // @todo here also variable usage dep. on top-level item?
}


MP2NL_Expr MP2NLModelAPI::AddExpression(const NLAffineExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_NLAffine); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const NLQuadExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_NLQuad); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AbsExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Abs); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const AndExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_And); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const OrExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Or); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const ExpExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Exp); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const LogExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Log); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const PowExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Pow); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const SinExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Sin); }
MP2NL_Expr MP2NLModelAPI::AddExpression(const CosExpression &expr)
{ return AddExpression(expr, ExpressionTypeID::ID_Cos); }


void MP2NLModelAPI::FinishProblemModificationPhase() { }




NLHeader MP2NLModelAPI::Header() {
  if (!hdr_is_current_) {
    PrepareModel();
    hdr_ = DoMakeHeader();
    hdr_is_current_ = true;
  }
  return hdr_;
}

void MP2NLModelAPI::PrepareModel() {
  MapExpressions();
  MarkVars();
  SortVars();
  MarkItems();
}

void MP2NLModelAPI::MapExpressions() {
  for (const auto& info: obj_info_)
    MapExprTreeFromItemInfo(info);
  for (const auto& info: alg_con_info_)
    MapExprTreeFromItemInfo(info);
  for (const auto& info: log_con_info_)
    MapExprTreeFromItemInfo(info);
  ResetIniExprRetrievedFlags();
}

void MP2NLModelAPI::MarkVars() {
  mark_data_.var_prior_.resize(var_lbs_.size());
  mark_data_.n_var_lin_bin_ = 0;
  mark_data_.n_var_lin_int_ = 0;
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
  hdr.num_logical_cons = log_con_info_.size();

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
  hdr.num_con_nonzeros = mark_data_.n_con_nz_;

  /** Number of nonzeros in all objective gradients. */
  hdr.num_obj_nonzeros = mark_data_.n_obj_nz_;

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

  hdr.format = storedOptions_.nl_format_text_
                   ? mp::NLHeader::TEXT : mp::NLHeader::BINARY;

  return hdr;
}

int MP2NLModelAPI::ObjType(int i) {
  switch (obj_info_[i].GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinearObjective:
  case StaticItemTypeID::ID_NLObjective:
    return obj::MIN==((LinearObjective*)(obj_info_[i].GetPItem()))->obj_sense()
               ? 0 : 1;
  default:
    MP_RAISE("Unknown objective type");
  }
}

template <class ObjGradWriterFactory>
void MP2NLModelAPI::FeedObjGradient(int i, ObjGradWriterFactory& svwf) {
  switch (obj_info_.at(i).GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinearObjective:
  case StaticItemTypeID::ID_NLObjective: {
    const auto& obj = *((LinearObjective*)(obj_info_[i].GetPItem()));
    if (obj.num_terms()) {
      auto svw = svwf.MakeVectorWriter(obj.num_terms());
      for (int j=0; j<obj.num_terms(); ++j)
        svw.Write(GetNewVarIndex( obj.GetLinTerms().var(j) ),      // new ordering
                  obj.GetLinTerms().coef(j));
    }
  } break;
  default:
    MP_RAISE("Unknown objective type");
  }
}

template <class ObjExprWriter>
void MP2NLModelAPI::FeedObjExpression(int iobj, ObjExprWriter& ew) {
  const auto& obj_info = obj_info_.at(iobj);
  switch (obj_info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinearObjective:
    ew.NPut(0.0);
    break;
  case StaticItemTypeID::ID_NLObjective: {
    const auto& obj = *((NLObjective*)(obj_info.GetPItem()));
    if (obj.HasExpr()) {
      FeedExpr(GetExpression(obj), ew);
    } else
      ew.NPut(0.0);
  } break;
  default:
    MP_RAISE("Unknown objective type");
  }
}

template <class VarBoundsWriter>
void MP2NLModelAPI::FeedVarBounds(VarBoundsWriter& vbw) {
  for (size_t i = 0; i < var_lbs_.size(); i++) {
    auto i_old = GetOldVarIndex(i);
    vbw.WriteLbUb(var_lbs_[i_old], var_ubs_[i_old]);
  }
}

template <class ConBoundsWriter>
void MP2NLModelAPI::FeedConBounds(ConBoundsWriter& cbw) {
  for (size_t i=0; i<alg_con_info_.size(); ++i) {          // no constraint permutations
    switch (alg_con_info_[i].GetStaticTypeID()) {
    case StaticItemTypeID::ID_LinConRange: {
      const auto& lcon = *((LinConRange*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_LinConLE: {
      const auto& lcon = *((LinConLE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_LinConEQ: {
      const auto& lcon = *((LinConEQ*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_LinConGE: {
      const auto& lcon = *((LinConGE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
    } break;
    case StaticItemTypeID::ID_NLConstraint: {
      const auto& lcon = *((NLConstraint*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{GetLower(lcon), GetUpper(lcon)} );
    } break;
    case StaticItemTypeID::ID_NLAssignLE: {  // -resvar + (expr) >= 0
      const auto& lcon = *((NLAssignLE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{0.0, Infinity()} );
    } break;
    case StaticItemTypeID::ID_NLAssignEQ: {
      const auto& lcon = *((NLAssignEQ*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{0.0, 0.0} );
    } break;
    case StaticItemTypeID::ID_NLAssignGE: {
      const auto& lcon = *((NLAssignGE*)(alg_con_info_[i].GetPItem()));
      cbw.WriteAlgConRange( AlgConRange{MinusInfinity(), 0.0} );
    } break;
    case StaticItemTypeID::ID_NLComplementarity: {
      const auto& lcon = *((NLComplementarity*)(alg_con_info_[i].GetPItem()));
      AlgConRange bnd;
      auto j = lcon.GetVariable();
      auto ct = lcon.GetExpression().constant_term(); // @todo NLExpression?
      bnd.L = MinusInfinity();
      bnd.U = Infinity();
      bnd.k = 0;
      if (var_lbs_[j] > MinusInfinity()) {
        bnd.k = 1;
        bnd.L = -ct;         // empty expr then
      }
      if (var_ubs_[j] < Infinity()) {
        bnd.k |= 2;
        bnd.U = -ct;
      }
      if (3==bnd.k) {
        bnd.L = MinusInfinity();
        bnd.U = Infinity();  // expr should be ct.
      }
      assert(bnd.k);
      bnd.cvar = GetNewVarIndex(j);
      cbw.WriteAlgConRange( bnd );
    } break;
    default:
      MP_RAISE("Unknown algebraic constraint type");
    }
  }
}

template <class ConLinearExprWriterFactory>
void MP2NLModelAPI::FeedLinearConExpr(int i, ConLinearExprWriterFactory& svwf) {
  const auto* pitem = alg_con_info_[i].GetPItem();
  switch (alg_con_info_[i].GetStaticTypeID()) {
  case StaticItemTypeID::ID_LinConRange: {
    FeedLinearConBody( *((LinConRange*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_LinConLE: {
    FeedLinearConBody( *((LinConLE*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_LinConEQ: {
    FeedLinearConBody( *((LinConEQ*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_LinConGE: {
    FeedLinearConBody( *((LinConGE*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_NLConstraint: {
    FeedLinearConBody(
        ((const NLConstraint*)(pitem))->GetMainCon(), svwf);
  } break;
  case StaticItemTypeID::ID_NLAssignLE: {
    FeedLinearConBody( *((NLAssignLE*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_NLAssignEQ: {
    FeedLinearConBody( *((NLAssignEQ*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_NLAssignGE: {
    FeedLinearConBody( *((NLAssignGE*)(pitem)), svwf);
  } break;
  case StaticItemTypeID::ID_NLComplementarity: {
    FeedLinearConBody(
        ((const NLComplementarity*)(pitem))->GetExpression(), svwf);
  } break;
  default:
    MP_RAISE("Unknown algebraic constraint type");
  }
}

template <class ExprBody, class ConLinearExprWriterFactory>
void MP2NLModelAPI::FeedLinearConBody(
    const ExprBody& algcon,
    ConLinearExprWriterFactory& svwf) {
  if (algcon.size()) {
    auto svw = svwf.MakeVectorWriter(algcon.size());
    for (int j=0; j<algcon.size(); ++j)
      svw.Write(GetNewVarIndex( algcon.var(j) ),      // new ordering
                algcon.coef(j));
  }
}

template <class ConExprWriter>
void MP2NLModelAPI::FeedConExpression(
    int icon, ConExprWriter& ew) {
  if (icon < (int)alg_con_info_.size())
    FeedAlgConExpression(icon, ew);
  else
    FeedLogicalConExpression(icon - alg_con_info_.size(), ew);
}

template <class ConExprWriter>
void MP2NLModelAPI::FeedAlgConExpression(
    int icon, ConExprWriter& ew) {
  // We could use a virtual function
  // to call GetExpression() instead of switch().
  const auto& con_info = alg_con_info_.at(icon);
  const auto* pitem = con_info.GetPItem();
  switch (con_info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_NLConstraint: {
    const auto& con = *((NLConstraint*)(pitem));
    if (con.HasExpr()) {
      FeedExpr(GetExpression(con), ew);
    } else
      ew.NPut(0.0);                   // linear case
  } break;
  case StaticItemTypeID::ID_NLAssignLE:
    FeedExpr(GetExpression(*((NLAssignLE*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLAssignEQ:
    FeedExpr(GetExpression(*((NLAssignEQ*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLAssignGE:
    FeedExpr(GetExpression(*((NLAssignGE*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLComplementarity: {
    const auto& cc = *((NLComplementarity*)(pitem));
    auto j = cc.GetVariable();
    if (var_lbs_[j]<=MinusInfinity()   // @todo proper expressions
        || var_ubs_[j]>=Infinity())
      ew.NPut(0.0);
    else
      ew.NPut(cc.GetExpression().constant_term());
  } break;
  default:
    ew.NPut(0.0);                     // must be linear constr
    break;
  case StaticItemTypeID::ID_NLLogical:
  case StaticItemTypeID::ID_NLImpl:
  case StaticItemTypeID::ID_NLEquivalence:
  case StaticItemTypeID::ID_NLRimpl:
  case StaticItemTypeID::ID_LinearObjective:
  case StaticItemTypeID::ID_NLObjective:
  case StaticItemTypeID::ID_None:
    MP_RAISE("Unknown algebraic constraint type");
  }
}

template <class ConExprWriter>
void MP2NLModelAPI::FeedLogicalConExpression(
    int icon, ConExprWriter& ew) {
  const auto& con_info = log_con_info_.at(icon);
  const auto* pitem = con_info.GetPItem();
  switch (con_info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_NLLogical:
    FeedExpr(GetExpression(*((NLLogical*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLImpl:
    FeedExpr(GetExpression(*((NLImpl*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLEquivalence:
    FeedExpr(GetExpression(*((NLEquivalence*)(pitem))), ew);
    break;
  case StaticItemTypeID::ID_NLRimpl:
    FeedExpr(GetExpression(*((NLRimpl*)(pitem))), ew);
    break;
  default:
    MP_RAISE("Unknown logical constraint type");
  }
}

template <class ExprWriter>
void MP2NLModelAPI::FeedExpr(Expr expr, ExprWriter& ew) {
  if (expr.IsEmptyExpr())
    ew.NPut(0.0);
  else if (expr.IsVariable()) {       // @todo def vars
    ew.VPut(expr.GetVarIndex(),
            GetVarName(expr.GetVarIndex()));
  } else
    FeedOpcode(expr, ew);
}

template <class ExprWriter>
void MP2NLModelAPI::FeedOpcode(Expr expr, ExprWriter& ew) {
  assert(expr.IsExpression());
  const auto& expr_info = expr_info_.at(expr.GetExprIndex());
  const auto* pitem = expr_info.GetPItem();
  switch (expr_info.GetExprTypeID()) {
  case ExpressionTypeID::ID_NLAffine:
    FeedAlgebraic(*(const NLAffineExpression*)pitem, ew);
    break;
  case ExpressionTypeID::ID_NLQuad:
    FeedAlgebraic(*(const NLQuadExpression*)pitem, ew);
    break;

  case ExpressionTypeID::ID_Abs:
    FeedOpcodeArgs(*(const AbsExpression*)pitem, ew.OPut1(nl::ABS));
    break;
  case ExpressionTypeID::ID_And:
    FeedOpcodeArgs(*(const AndExpression*)pitem,
                   ew.OPutN(nl::AND,
                            GetNumArguments(*(const AndExpression*)pitem)));
    break;
  case ExpressionTypeID::ID_Or:
    FeedOpcodeArgs(*(const OrExpression*)pitem,
                   ew.OPutN(nl::OR,
                            GetNumArguments(*(const OrExpression*)pitem)));
    break;
  case ExpressionTypeID::ID_Exp:
    FeedOpcodeArgs(*(const ExpExpression*)pitem, ew.OPut1(nl::EXP));
    break;
  case ExpressionTypeID::ID_Log:
    FeedOpcodeArgs(*(const LogExpression*)pitem, ew.OPut1(nl::LOG));
    break;
  case ExpressionTypeID::ID_Pow:
    FeedOpcodeArgs(*(const PowExpression*)pitem, ew.OPut2(nl::POW_CONST_EXP));
    break;
  case ExpressionTypeID::ID_Sin:
    FeedOpcodeArgs(*(const SinExpression*)pitem, ew.OPut1(nl::SIN));
    break;
  case ExpressionTypeID::ID_Cos:
    FeedOpcodeArgs(*(const CosExpression*)pitem, ew.OPut1(nl::COS));
    break;
  default:
    MP_RAISE("MP2NL: unknown expression type");
  }
}

template <class AlgMPExpr, class ExprWriter>
void MP2NLModelAPI::FeedAlgebraic(
    const AlgMPExpr& e, ExprWriter& ew) {
  int n_args = GetLinSize(e) + GetQuadSize(e) + bool(GetConstTerm(e));
  auto write_args = [&,this](const auto& e, auto ew_args0) {
    for (int i=0; i<GetLinSize(e); ++i) {
      if (1.0==GetLinCoef(e, i))
        ew_args0.EPut(GetLinTerm(e, i));
      else {
        auto ew_args1 = ew_args0.OPut2(nl::MUL);
        ew_args1.NPut(GetLinCoef(e, i));
        ew_args1.EPut(GetLinTerm(e, i));
      }
    }
    if constexpr (std::is_same_v<AlgMPExpr, NLQuadExpression>) {
      for (int i=0; i<GetQuadSize(e); ++i) {
        if (1.0==GetQuadCoef(e, i)) {
          auto ew_args1 = ew_args0.OPut2(nl::MUL);
          ew_args1.EPut(GetQuadTerm1(e, i));
          ew_args1.EPut(GetQuadTerm2(e, i));
        } else {
          auto ew_args1 = ew_args0.OPut2(nl::MUL);
          ew_args1.NPut(GetQuadCoef(e, i));
          auto ew_args2 = ew_args1.OPut2(nl::MUL);
          ew_args2.EPut(GetQuadTerm1(e, i));
          ew_args2.EPut(GetQuadTerm2(e, i));
        }
      }
    }
    if (GetConstTerm(e))
      ew_args0.NPut(GetConstTerm(e));
  };
  if (!n_args)
    ew.NPut(0.0);
  else if (1==n_args)
    write_args(e, std::move(ew));
  else if (2==n_args)
    write_args(e, ew.OPut2(nl::ADD));
  else
    write_args(e, ew.OPutN(nl::SUM, n_args));
}

template <class MPExpr, class ArgWriter>
void MP2NLModelAPI::FeedOpcodeArgs(
    const MPExpr& e, ArgWriter ew_arg) {
  for (int i=0; i<GetNumArguments(e); ++i)
    ew_arg.EPut(GetArgExpression(e, i));
  for (int i=0; i<GetNumParameters(e); ++i)
    ew_arg.NPut(GetParameter(e, i));
}


template <class ColSizeWriter>
void MP2NLModelAPI::FeedColumnSizes(ColSizeWriter& csw) {
  if (WantColumnSizes())
    for (int i=0; i < var_lbs_.size()-1; ++i)        // use old ordering
      csw.Write(mark_data_.col_sizes_orig_[ GetOldVarIndex( i ) ]);
}


/** Initial primal guesses.
   *
   *  Implementation: write all meaningfuls entries (incl. zeros.)
   *      if (ini_guess.size()) {
   *        auto ig = igw.MakeVectorWriter(ini_guess.size());
   *        for (size_t i=0; i<ini_guess.size(); ++i)
   *          ig.Write(ini_guess[i].index_, ini_guess[i].value_);
   *      }
   */
template <class IGWriter>
void MP2NLModelAPI::FeedInitialGuesses(IGWriter& igw) {
  auto x0 = GetNLSolver().GetCallbacks()->GetInitialGuesses();
  if (x0.size()) {
    auto ig = igw.MakeVectorWriter(x0.size());
    for (size_t i=0; i<x0.size(); ++i) {
      ig.Write(GetNewVarIndex(x0[i].first), x0[i].second);
    }
  }
}

/** Initial dual guesses. */
template <class IDGWriter>
void MP2NLModelAPI::FeedInitialDualGuesses(IDGWriter& igw) {
  auto y0 = GetNLSolver().GetCallbacks()->GetInitialDualGuesses();
  if (y0.size()) {
    auto ig = igw.MakeVectorWriter(y0.size());
    for (size_t i=0; i<y0.size(); ++i) {
      ig.Write(i, y0[i]);
    }
  }
}


/** Feed suffixes.
     *
     *  For constraints, assume ordering:
     *  first algebraic, then logical.
   *
   *  Implementation: write all non-0 entries (0 is the default.)
   *      while (....) {
   *        auto sw = swf.StartIntSuffix(  // or ...DblSuffix
   *          suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          sw.Write(index[i], value[i]);
   *      }
     */
template <class SuffixWriterFactory>
void MP2NLModelAPI::FeedSuffixes(SuffixWriterFactory& swf) {
  // @todo SOS constraints
  auto p_qc = p_nlsi_->GetCallbacks();
  auto suffixnames = p_qc->GetSuffixNames();
  auto write1suf = [this](auto& sw, int kind, const auto& vals) {
    for (size_t i=0; i<vals.size(); ++i) {
      if (auto val = vals[i]) {
        int i0=i;
        if (suf::VAR==kind)
          i0 = this->GetNewVarIndex(i);
        sw.Write(i0, val);
      }
    }
  };
  for (const auto& sufname: suffixnames) {
    auto modelsuffix = p_qc->GetModelSuffix(sufname);
    for (int kind=0; kind<modelsuffix.values_.size(); ++kind) {
      const auto& vals = modelsuffix.values_[kind];
      if (vals.size()) {                 // even all-0 suffixes
        auto nnz = std::count_if(vals.begin(), vals.end(),
                                 [](auto n){ return bool(n); });
        if (modelsuffix.flags_ & suf::FLOAT) {
          auto sw = swf.StartDblSuffix(sufname.c_str(),
                                       kind | suf::FLOAT, nnz);
          write1suf(sw, kind, vals);
        } else {
          auto sw = swf.StartIntSuffix(sufname.c_str(), kind, nnz);
          write1suf(sw, kind, vals);
        }
      }
    }
  }
}

template <class ColNameWriter>
void MP2NLModelAPI::FeedRowAndObjNames(ColNameWriter& wrt) {
  auto has_name0 = [this](const auto& infos) {
    return infos.size() && infos.front().GetDispatcher().GetName(infos.front().GetPItem());
  };
  if ((has_name0(alg_con_info_)
       || has_name0(log_con_info_) || has_name0(obj_info_)) && wrt) {
    auto write_names = [this,&wrt](const auto& infos) {
      for (size_t i=0; i<infos.size(); ++i) {
        const auto* nm
            = infos[i].GetDispatcher().GetName(infos[i].GetPItem());
        wrt << (nm ? nm : "..");
      }
    };
    write_names(alg_con_info_);    // @todo any permutations
    write_names(log_con_info_);
    write_names(obj_info_);
  }
}

template <class ColNameWriter>
void MP2NLModelAPI::FeedColNames(ColNameWriter& wrt) {
  if (var_names_.size() && wrt) {
    assert(var_names_.size() == var_lbs_.size());
    for (size_t i=0; i<var_names_.size(); ++i)
      wrt << var_names_[ GetOldVarIndex(i) ];
  }
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
    nlsol_.SetFileStub(mapi_.GetFileStub());
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

  /// Suffix names
  std::set<std::string> GetSuffixNames() override {
    std::set<std::string> result;
    for (const auto& suf: suf_map_)
      result.insert(suf.first);
    return result;
  }

  /// Get model suffix with given name
  const MP2NLModelSuffix& GetModelSuffix(
      const std::string& name) override {
    return suf_map_.at(name);
  }

  /// Num alg cons
  int GetNumAlgCons() const override {
    return n_alg_con_;
  }

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
  int OnAMPLOptions(const AMPLOptions& ) {
    suf_map_.clear();              // clear for new solution
    return 0;
  }

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
    if (int nac_sol = rd.Size()) {
      auto n_alg_cons = Header().num_algebraic_cons;
      if (nac_sol > n_alg_cons) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_more_alg_cons",
            "The subsolver reported more duals than algebraic constraints");
      } else if (nac_sol < n_alg_cons) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_fewer_alg_cons",
            "The subsolver reported fewer duals than algebraic constraints");
      }
      duals_.reserve(nac_sol);
      while (rd.Size())
        duals_.push_back( rd.ReadNext() );                  // no cperm
    }
  }

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  template <class VecReader>
  void OnPrimalSolution(VecReader& rd) {
    primals_.clear();
    if (int nv_sol = rd.Size()) {
      auto n_vars = Header().num_vars;
      if (nv_sol > n_vars) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_more_vars",
            "The subsolver reported more variables");
      } else if (nv_sol < n_vars) {
        mapi_.GetEnv().AddWarning(
            "MP2NL_subsolver_solution_fewer_vars",
            "The subsolver reported fewer variables");
      }
      primals_.resize(n_vars);
      int j=0;
      for ( ; rd.Size(); ++j ) {
        int j0 = mapi_.GetOldVarIndex(j);
        assert(j0>=0 && j0 < n_vars);
        primals_[j0] = rd.ReadNext();
      }
    }
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
    OnSuffix(sr);
  }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
  template <class SuffixReader>
  void OnDblSuffix(SuffixReader& sr) {
    OnSuffix(sr);
  }


protected:
  template <class SuffixReader>
  void OnSuffix(SuffixReader& sr) {
    if (!if_suf_data_registered_) {
      if_suf_data_registered_ = true;
      RegisterSuffixData();
    }
    const auto& si = sr.SufInfo();
    int kind = si.Kind();
    int nmax = nitems_[kind & 3];
    const std::string& name = si.Name();
    const std::string& table = si.Table();
    auto& modelsuf = suf_map_[name];
    modelsuf.name_ = name;
    // set to FLOAT if at least one
    modelsuf.flags_ |= (kind & suf::FLOAT);
    if (modelsuf.table_.size() < table.size())
      modelsuf.table_ = table;
    auto& suf = modelsuf.values_.at(kind & 3);
    suf.clear();
    suf.resize(nitems_[kind & 3]);
    while (sr.Size()) {
      auto sparse_entry = sr.ReadNext();
      if (sparse_entry.first<0 || sparse_entry.first>=nmax) {
        sr.SetError(NLW2_SOLRead_Bad_Suffix,
                    "bad suffix element index");
        return;
      }
      int i0 = sparse_entry.first;
      if (0 == (kind & 4))                   // variable suffix
        i0 = mapi_.GetOldVarIndex(i0);
      suf[i0] = sparse_entry.second;
    }
  }

  /// Register some data for suffix reporting
  void RegisterSuffixData() {
    const auto& hdr = Header();

    nitems_[0] = hdr.num_vars;
    nitems_[1] = hdr.num_algebraic_cons + hdr.num_logical_cons;
    nitems_[2] = hdr.num_objs;
    nitems_[3] = 1;            // N problems

    n_alg_con_ = hdr.num_algebraic_cons;
  }


private:
  MP2NLModelAPI& mapi_;

  mp::NLUtils utils_;

  mp::NLSolver nlsol_ { &utils_ };

  bool if_suf_data_registered_ {};
  std::array<int, 4> nitems_ {};
  int n_alg_con_ {};               // For Backend to split alg + log cons

  std::unordered_map<std::string, MP2NLModelSuffix> suf_map_;

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

template <class Lambda>
void MP2NLModelAPI::DispatchStaticItem(
    const ItemInfo& info, Lambda lambda) {
  switch (info.GetStaticTypeID()) {
  case StaticItemTypeID::ID_None:
    MP_RAISE("MP2NL static item type not specified");
  case StaticItemTypeID::ID_LinearObjective:
    lambda(*(const LinearObjective*)info.GetPItem());
    break;
  default:
    MP_RAISE("MP2NL: unknown static item type");
  }
        /*
      ID_NLObjective,

      ID_LinConRange,
      ID_LinConLE,
      ID_LinConEQ,
      ID_LinConGE,
      ID_NLConstraint,
      ID_NLAssignLE,
      ID_NLAssignEQ,
      ID_NLAssignGE,
      ID_NLLogical,
      ID_NLEquivalence,
      ID_NLImpl,
      ID_NLRimpl,
      ID_IndicatorConstraintLinLE,
      ID_IndicatorConstraintLinEQ,
      ID_IndicatorConstraintLinGE,
      ID_SOS1Constraint,
      ID_SOS2Constraint,
      ID_NLComplementarity  */
}

void MP2NLModelAPI::CreateInterfaces() {
  p_nls_ = std::make_unique<MP2NLSolverImpl>(*this);
  this->MP2NLCommonInfo::p_nlsi_ = p_nls_.get();
}

} // namespace mp
