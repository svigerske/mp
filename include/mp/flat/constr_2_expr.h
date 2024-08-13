#ifndef CONSTR_2_EXPR_H
#define CONSTR_2_EXPR_H

#include <functional>
#include <type_traits>

#include "mp/flat/constr_keeper.h"
#include "mp/flat/nl_expr/constr_nl.h"
#include "mp/valcvt-link.h"

namespace mp {

/// A mix-in base class
/// facilitating "inlining" of some functional constraints
/// into expression trees.
template <class Impl>
class Constraints2Expr {
public:
  /// Mark some functional constraints as expressions.
  /// Then convert some objectives and constraints
  /// to "NL...." items using the expressions.
  /// The auxiliary result variables of the functional
  /// constraints which were marked as expressions,
  /// will be handled in a special way.
  void Convert2NL() {
    MPD( MarkExpressions() );
    MPD( GetModel() ).ConvertAllWithExpressions(*(Impl*)this);
  }

  /// Mark which functional constraints to be used as expressions,
  /// vs assigning their result to a variable
  void MarkExpressions() {
    MPD( MarkAllResultVarsAsVars() );      // tentatively
    MPD( GetModel() ).MarkExprResultVars(*(Impl*)this);
    MPD( GetModel() ).MarkArguments(*(Impl*)this);
  }

  /// Consider marking the result variables as
  /// possible expressions
  template <class Con>
  void ConsiderMarkingResultVar(
      const Con& con, int , ExpressionAcceptanceLevel eal) {
    assert(ExpressionAcceptanceLevel::NotAccepted!=eal);
    if (con.HasResultVar()) {         // A functional constraint
      assert(                     // Check: the result var has \a con as the init expr
          MPD( template GetInitExpressionOfType<Con>(con.GetResultVar()) )
          == &con);
      MPD( MarkAsExpression(con.GetResultVar()) );   // can be changed later
    }
  }

  /// Consider marking the argument variables as
  /// "explicit variables" (not expressions.)
  /// Generic request.
  /// For func cons, once not accepted as expressions, consider their args
  /// as variables.
  /// For static cons, always.
  template <class Con>
  void ConsiderMarkingArguments(
      const Con& con, int i, ExpressionAcceptanceLevel eal) {
    bool fMarkArgs = false;
    if (con.HasResultVar())    // func cons: non-accepted ones by default
      fMarkArgs = (ExpressionAcceptanceLevel::NotAccepted==eal);
    else
      fMarkArgs = true;        // static cons: all non-algebraic by default
    if (fMarkArgs)
      MPD( DoMarkArgsAsVars(con, i) );
  }


  /// Convert algebraic constraint for use with expressions
  /// @param i: constraint index
  /// @param cal: (user chosen) acceptance level for \a con.
  ///   While it is recommended that flat algebraic constraints
  ///   are still accepted along with NLConstraint,
  ///   user might set acc:linle=0 etc.
  /// @return true if this constraint has been eliminated/replaced.
  template <class Body, class RhsOrRange>
  bool ConvertWithExpressions(
      const AlgebraicConstraint<Body, RhsOrRange>& con,
      int i,
      ConstraintAcceptanceLevel cal) {
    /// Replace \a con by a NLConstraint,
    /// if either the ModelAPI does not accept it,
    /// or the linear/quadratic terms have expressions
    /// (and then they are non-flat.)
    if (ConstraintAcceptanceLevel::Recommended != cal
        || HasExpressionArgs(con.GetBody())) {
      ConvertToNLCon(con, i);
      return true;                              // to remove the original \a con
    }
    return false;
  }

  /// Convert conditional algebraic constraint for use with expressions.
  /// Basically we need the LHS to be an expression.
  /// This can be simplified if we switch to do it earlier in flat conversions.
  /// @return true if this constraint has been eliminated/replaced.
  template <class Body, class RhsOrRange>
  bool ConvertWithExpressions(
      const ConditionalConstraint< AlgebraicConstraint<Body, RhsOrRange> >& con,
      int i,
      ConstraintAcceptanceLevel ) {
    ConsiderExplicifyingExpression(con, i);        // this is a func con too
    if (!con.GetConstraint().GetBody().is_variable()) {      // already a variable
      ConvertConditionalConLHS(con, i);
      return true;
    }
    return false;
  }

  /// Convert other functional constraints for use with expressions.
  /// Basically, see if the result variable is explicit.
  /// @return true if this constraint has been eliminated/replaced.
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<FunctionalConstraint, FuncCon>, bool > = true >
  bool ConvertWithExpressions(
      const FuncCon& con, int i, ConstraintAcceptanceLevel ) {
    ConsiderExplicifyingExpression(con, i);
    return false;                          // leave it active
  }

  /// Convert complementarity constraint for use with expressions.
  /// Similarly to Conditional, we need the expression part to be an NL expression.
  /// @return true if this constraint has been eliminated/replaced.
  template <class Expr>
  bool ConvertWithExpressions(
      const ComplementarityConstraint<Expr>& con,
      int i,
      ConstraintAcceptanceLevel ) {
    if (!con.GetExpression().is_variable()) {             // already a variable
      ConvertComplementarityExpr(con, i);
      return true;
    }
    return false;
  }

  /// Template for Indicators: do nothing.
  /// But check that they are flat?
  template <class SubCon>
  bool ConvertWithExpressions(
      const IndicatorConstraint<SubCon>& , int , ConstraintAcceptanceLevel ) {
    return false;
  }

  /// SOS1: do nothing.
  /// But check that they are flat?
  bool ConvertWithExpressions(
      const SOS1Constraint& con, int , ConstraintAcceptanceLevel ) {
    for (int v: con.GetArguments()) {
      assert(MPCD( IsProperVar(v) ));
    }
    return false;
  }

  /// SOS1: do nothing.
  /// But check that they are flat?
  bool ConvertWithExpressions(
      const SOS2Constraint& con, int , ConstraintAcceptanceLevel ) {
    for (int v: con.GetArguments()) {
      assert(MPCD( IsProperVar(v) ));
    }
    return false;
  }

  /// NLConstraint: just produced.
  bool ConvertWithExpressions(
      const NLConstraint& , int , ConstraintAcceptanceLevel ) {
    return false;
  }
  /// NLLogical: just produced.
  bool ConvertWithExpressions(
      const NLLogical& , int , ConstraintAcceptanceLevel ) {
    return false;
  }
  /// NLEquivalence: just produced.
  bool ConvertWithExpressions(
      const NLEquivalence& , int , ConstraintAcceptanceLevel ) {
    return false;
  }
  /// NLImpl: just produced.
  bool ConvertWithExpressions(
      const NLImpl& , int , ConstraintAcceptanceLevel ) {
    return false;
  }
  /// NLRimpl: just produced.
  bool ConvertWithExpressions(
      const NLRimpl& , int , ConstraintAcceptanceLevel ) {
    return false;
  }


protected:
  /// Algebraic cons: no marking (when NLConstraint accepted?)
  /// @todo Do we need to consider NLConstraint / NLObjective at this step?
  /// Are they added during marking?
  template <class Body, class RhsOrRange>
  void DoMarkArgsAsVars(   // needs to appear before the most generic template
      const AlgebraicConstraint<Body, RhsOrRange>& , int ) { }

  /// Complementarity: only mark the var
  /// (NL accepts expressions for the expression part)
  template <class Expr>
  void DoMarkArgsAsVars(   // needs to appear before the most generic template
      const ComplementarityConstraint<Expr>& cc, int ) {
    MPD( MarkAsResultVar(cc.GetVariable()) );
  }

  /// Generic arguments marking call
  template <class Con>
  void DoMarkArgsAsVars(const Con& con, int ) {
    VisitArguments(con, MarkVar_);
  }


  /// Check if the algebraic body has expressions
  bool HasExpressionArgs(const QuadAndLinTerms& qlt) const {
    return HasExpressionArgs(qlt.GetLinTerms())
           || HasExpressionArgs(qlt.GetQPTerms());
  }

  bool HasExpressionArgs(const LinTerms& lt) const {
    for (auto v: lt.vars())
      if (!MPCD( IsProperVar(v) ))
        return true;
    return false;
  }

  bool HasExpressionArgs(const QuadTerms& qt) const {
    for (auto v: qt.vars1())
      if (!MPCD( IsProperVar(v) ))
        return true;
    for (auto v: qt.vars2())
      if (!MPCD( IsProperVar(v) ))
        return true;
    return false;
  }

  /// Convert algebraic con to a NLConstraint
  template <class Body, class RhsOrRange>
  void ConvertToNLCon(
      const AlgebraicConstraint<Body, RhsOrRange>& con, int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );   // link from \a con
    LinTerms lt;
    /// exprTerm will be a LinearFunctionalConstraint or a Quadratic...
    auto exprTerm = ExtractLinAndExprArgs(con.GetBody(), lt);
    assert(0.0 == exprTerm.GetArguments().constant_term());
    /// Store full LFC only if it is not 1.0*var
    int exprResVar = -1;
    if (exprTerm.GetArguments().is_variable()) {
      exprResVar = exprTerm.GetArguments().get_representing_variable();
      assert( !MPCD( IsProperVar(exprResVar) ) );     // is an expr
    } else {
      assert( !exprTerm.GetArguments().empty() );     // has more terms, or coef != 1.0
      exprResVar = MPD( AssignResultVar2Args(std::move(exprTerm)) );
      if (!MPCD(VarHasMarking(exprResVar)))           // mark as expr if new
        MPD( MarkAsExpression(exprResVar) );
      /// Exists and marked a variable
      else if (MPCD( IsProperVar(exprResVar) )) {     // Not an expression after all
        lt.add_term(1.0, exprResVar);
        exprResVar = -1;                              // no expression
        lt.sort_terms();
      }
    }  // else, stays as -1
    NLConstraint nlc{lt, exprResVar,
                     { con.GetRhsOrRange().lb(), con.GetRhsOrRange().ub() },
                     false};                  // no sorting
    MPD( AddConstraint( std::move(nlc) ) );
  }

  /// Extract linear and expression args
  /// from Quad+Linear terms
  /// @param ltout: pure linear part (coefs*vars)
  ///   for NLConstraint
  /// @return expression part (QFC)
  QuadraticFunctionalConstraint ExtractLinAndExprArgs(
      const QuadAndLinTerms& qltin,
      LinTerms& ltout) {
    /// Create full QFC because we have QP terms
    LinTerms lt_in_expr = SplitLinTerms(qltin.GetLinTerms(), ltout);
    assert(qltin.GetQPTerms().size());
    return
        {{{std::move(lt_in_expr), qltin.GetQPTerms()}, 0.0}};
  }

  /// Extract linear and expression args
  /// from Linear terms.
  /// @param ltout: pure linear part (coefs*vars)
  ///   for NLConstraint
  /// @return expression part (LFC)
  LinearFunctionalConstraint ExtractLinAndExprArgs(
      const LinTerms& ltin,
      LinTerms& ltout) {
    LinTerms lt_in_expr = SplitLinTerms(ltin.GetLinTerms(), ltout);
    return {{std::move(lt_in_expr), 0.0}};
  }

  /// Split linterms into pure-var and expressions
  /// @return the expressions part
  LinTerms SplitLinTerms(const LinTerms& ltin, LinTerms lt_out_vars) {
    LinTerms result;
    int nvars=0;
    for (int v: ltin.vars())
      if (MPCD( IsProperVar(v) ))
        ++nvars;
    lt_out_vars.reserve(nvars);
    result.reserve(ltin.size() - nvars);
    int v=0;
    for (size_t i=0; i<ltin.size(); ++i) {
      if (MPCD( IsProperVar(v = ltin.var(i)) ))
        lt_out_vars.add_term(ltin.coef(i), v);
      else
        result.add_term(ltin.coef(i), v);
    }
    return result;
  }

  /// Convert the LHS of a conditional con
  template <class Body, class RhsOrRange>
  void ConvertConditionalConLHS(
      const ConditionalConstraint< AlgebraicConstraint<Body, RhsOrRange> >& con,
      int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    /// Create a functional constraint from the LHS
    auto fc = MakeFunctionalConstraint(
        AlgebraicExpression{con.GetConstraint().GetArguments(), 0.0});
    auto resvar = MPD( AssignResultVar2Args(std::move(fc)) );
    /// resvar can be a proper variable - ModelAPI should flexibly handle this
    LinTerms lt { {1.0}, {resvar} };
    ConditionalConstraint< AlgebraicConstraint<LinTerms, RhsOrRange> >
        ccnew { { std::move(lt), con.GetConstraint().GetRhsOrRange() } };
    MPD( RedefineVariable(con.GetResultVar(), std::move(ccnew)) );  // Use new CondCon
  }

  /// Convert the expression part of complementarity.
  /// Similar to the argument of a conditional con.
  template <class Expr>
  void ConvertComplementarityExpr(
      const ComplementarityConstraint<Expr>& con,
      int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    /// Create a functional constraint from the LHS
    auto fc = MakeFunctionalConstraint(con.GetExpression());
    auto resvar = MPD( AssignResultVar2Args(std::move(fc)) );
    /// resvar can be a proper variable - ModelAPI should flexibly handle this
    LinTerms lt { {1.0}, {resvar} };
    ComplementarityConstraint< AlgebraicExpression<LinTerms> >
        ccnew { AffineExpr{ std::move(lt), 0.0 }, con.GetVariable() };
    MPD( AddConstraint(std::move(ccnew)) );  // Use new CondCon
  }

  /// Consider explicifying an expression
  template <class FuncCon>
  void ConsiderExplicifyingExpression(const FuncCon& con, int i) {
    if (MPCD( IsProperVar(con.GetResultVar()) ))
      DoExplicify(con, i);
  }

  /// Add expr = var assignment for algebraic expression
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<NumericFunctionalConstraintTraits, FuncCon>, bool > = true >
  void DoExplicify(const FuncCon& con, int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    assert(!con.GetContext().IsNone());
    auto resvar = con.GetResultVar();
    AlgConRange rng {-INFINITY, 0.0};                     // ctx-: var >= expr(var)
    if (con.GetContext().IsMixed())
      rng = {0.0, 0.0};
    else if (con.GetContext().HasPositive())
      rng = {0.0, INFINITY};
    MPD( AddConstraint(                                   // -var + expr (in) rng
        NLConstraint{ { {-1.0}, {resvar} }, resvar, rng } ) );
  }

  /// Add expr = var assignment for logical expression (NLEquivalence, NLImpl, NLRImpl).
  /// If this constraint is root-level (is true),
  /// add root explicifier instead (NLLogical).
  template <class FuncCon,
           std::enable_if_t<
               std::is_base_of_v<LogicalFunctionalConstraintTraits, FuncCon>, bool > = true >
  void DoExplicify(const FuncCon& con, int i) {
    auto alscope = MPD( MakeAutoLinker( con, i ) );       // link from \a con
    auto resvar = con.GetResultVar();
    if (MPCD( is_fixed(resvar) )                          // "root" context
        && MPCD( fixed_value(resvar) )) {                 // static logical constraint
      MPD( AddConstraint(NLLogical(resvar)) );
    }
    else {                                                // A (half-)reified function constraint
      assert(!con.GetContext().IsNone());
      if (con.GetContext().IsMixed())
        MPD( AddConstraint(NLEquivalence(resvar)) );
      else if (con.GetContext().HasPositive())
        MPD( AddConstraint(NLImpl(resvar)) );
      else
        MPD( AddConstraint(NLRimpl(resvar)) );
    }
  }


private:
  /// (Argument) variable visitor
  std::function<void( int )> MarkVar_ = [this](int v){
    MPD( MarkAsResultVar(v) );       // no recursion
  };
};

}  // namespace mp

#endif // CONSTR_2_EXPR_H
