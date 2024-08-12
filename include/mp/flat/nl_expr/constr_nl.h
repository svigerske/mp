#ifndef CONSTR_NL_H
#define CONSTR_NL_H

#include <cassert>

#include "mp/flat/constr_std.h"

namespace mp {

/// Class NLConstraint.
/// Algebraic range constraint with a linear part
/// and an expression term.
class NLConstraint
    : public LinConRange {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    return "NLConstraint";
  }

  /// Constructor.
  /// @param linexpr: linear part
  /// @param expr: result variable of the expression part
  /// @param rng: {lb, ub}
  /// @param fSort=true: whether to sort linear terms
  NLConstraint(
      const LinTerms& lt, int expr, AlgConRange rng,
      bool fSort=true)
      : LinConRange(lt, rng, fSort), expr_(expr) { }

  /// Has expression term?
  bool HasExpr() const { return expr_>=0; }

  /// Expression index.
  /// @note ModelAPI should call self.HasExpression
  ///   and self.GetExpression()
  ///   to obtain the expression term.
  int ExprIndex() const { assert(HasExpr()); return expr_; }

  /// Compute violation.
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x, bool logical=false) const {
    double bd = GetBody().ComputeValue(x);
    if (HasExpr())
      bd += x[ExprIndex()];     // Add expr value. Assume it's precomputed
    if (!logical) {
      if (lb() > bd)
        return {lb() - bd, lb()};
      if (bd > ub())
        return {bd - ub(), ub()};
      return
          {std::max( // negative. Same for strict cmp?
               lb() - bd, bd - ub()),
           0.0};
    }
    return {double(!is_valid(bd)), 1.0};
  }


private:
  int expr_ {-1};
};


/// Export to JSON
inline void WriteJSON(JSONW jw,
                      const NLConstraint& nlc) {
  WriteJSON(jw["lin_part"], nlc.GetBody());
  if (nlc.HasExpr())
    WriteJSON(jw["expr_index"], nlc.ExprIndex());
  WriteJSON(jw["rhs_or_range"], nlc.GetRhsOrRange());
}

}  // namespace mp

#endif // CONSTR_NL_H
