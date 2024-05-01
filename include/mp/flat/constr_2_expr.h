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

  /// Consider marking the result and argument variables as
  /// "explicit variables" (not expressions)
  template <class Con>
  void ConsiderMarkingResultAndArgVars(
      const Con& con, int i, ExpressionAcceptanceLevel eal) {
    if (con.HasResultVar()) {         // A functional constraint
      if (ExpressionAcceptanceLevel::NotAccepted==eal) {
        MPD( MarkAsResultVar(con.GetResultVar()) );
      }
    }                          // Any constraint
    MPD( ConsiderMarkingArgumentsAsVars(con, i, eal) );
  }

  /// Generic request to consider marking arguments
  template <class Con>
  void ConsiderMarkingArgumentsAsVars(
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
