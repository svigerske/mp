#ifndef CONSTR_2_EXPR_H
#define CONSTR_2_EXPR_H

#include <functional>

#include "mp/flat/constr_keeper.h"

namespace mp {

/// A mix-in base class
/// facilitating "inlining" of some functional constraints
/// into expression trees.
template <class Impl>
class Constraints2Expr {
public:
  /// Convert some functional constraints to expressions
  void Convert2NL() {
    MPD( MarkExpressions() );
    MPD( GetModel() ).ConvertWithExpressions(*(Impl*)this);
  }

  /// Mark which functional constraints to be used as expressions,
  /// vs assigning their result to a variable
  void MarkExpressions() {
    MPD( MarkAllResultVarsAsVars() );
    MPD( GetModel() ).MarkExprResultVars(*(Impl*)this);
    MPD( GetModel() ).MarkArguments(*(Impl*)this);
  }

  /// Consider marking the result variables as
  /// possible expressions
  template <class Con>
  void ConsiderMarkingResultVar(
      const Con& con, int , ExpressionAcceptanceLevel eal) {
    if (con.HasResultVar()) {         // A functional constraint
      // Check that it will be expression, but possibly with a dedicated result variable
      if (ExpressionAcceptanceLevel::NotAccepted!=eal) {
        assert(                     // Check: the result var has \a con as the init expr
            MPD( template GetInitExpressionOfType<Con>(con.GetResultVar()) )
                       == &con);
        MPD( MarkAsExpression(con.GetResultVar()) );   // can be changed later
      }
    }
  }

  /// Consider marking the argument variables as
  /// "explicit variables" (not expressions.)
  /// Generic request.
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


protected:
  /// Algebraic cons: no marking (when NLConstraint accepted?)
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


private:
  /// (Argument) variable visitor
  std::function<void( int )> MarkVar_ = [this](int v){
    MPD( MarkAsResultVar(v) );       // no recursion
  };
};

}  // namespace mp

#endif // CONSTR_2_EXPR_H
