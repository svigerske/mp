/*
 Basic expression-based model API definitions.

 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov <gleb@ampl.com>
*/
#ifndef MODEL_API_BASE_H
#define MODEL_API_BASE_H

#include <utility>

#include "mp/flat/model_api_base.h"
#include "mp/flat/nl_expr/constr_nl.h"

namespace mp {

/// ModelAPIs handling expression trees should derive from
/// BasicExprModelAPI<Impl, Expr>,
/// with Impl the final implementation class (CRTP)
/// and Expr default-constructible.
template <class Impl, class ExprType=void*>
class BasicExprModelAPI
    : public BasicFlatModelAPI {
public:
  using Expr = ExprType;
  /// Placeholder for GetTypeName()
  static const char* GetTypeName()    { return "BasicExprModelAPI"; }

/// A ModelAPI accepting NL trees should declare this.
///
/// - NotAccepted: not compiled
/// - AcceptedButNotRecommended: compiled but off by default (option acc:_expr)
/// - Recommended: on by default
#define ACCEPT_EXPRESSION_INTERFACE(val) \
  static constexpr ExpressionAcceptanceLevel \
  ExpressionInterfaceAcceptanceLevel() { return ExpressionAcceptanceLevel::val; }

  /// Whether accepts NLObjective.
  /// SCIP does not.
  static int AcceptsNLObj() { return 1; }

/// Reuse inherited names
  USE_BASE_CONSTRAINT_HANDLERS(BasicFlatModelAPI)


  /// Placeholder for AddExpression<>()
  template <class Expression>
  Expr AddExpression(const Expression& e) {
    MP_RAISE(
        std::string("Not handling expression type '") +
        e.GetTypeName() +
        "'. Provide a handler or a converter method");
  }

/// Derived backends have to tell C++ to use default handlers if they are needed
/// when they overload AddExpression(), due to C++ name hiding
#define USE_BASE_EXPRESSION_HANDLERS(BaseBackend) \
  using BaseBackend::Expr; \
  using BaseBackend::AddExpression;


  /// NL model item accessors

  /// Get num linear terms
  int GetLinSize(const NLConstraint& nlc) const {
    return nlc.GetMainCon().size();
  }

  /// Get linear coef \a i.
  double GetLinCoef(const NLConstraint& nlc, int i) const {
    return nlc.GetMainCon().coef(i);
  }

  /// Get linear part var \a i.
  int GetLinVar(const NLConstraint& nlc, int i) const {
    assert(IsVarProper(i));
    return nlc.GetMainCon().var(i);
  }

  /// Get the expression term of an \a NLConstraint.
  /// @note Can return the dummy expression
  ///   via Impl::GetZeroExpression().
  ExprType GetExpression(const NLConstraint& nlc) {
    const auto i_expr = nlc.HasExpr()
                            ? nlc.ExprIndex() : -1;
    if (i_expr<0)
      return MPD( GetZeroExpression() );
    return GetPureInitExpression(i_expr);
  }

  /// Get NLConstraint's lower bound
  double GetLower(const NLConstraint& nlc) const {
    return nlc.GetMainCon().lb();
  }

  /// Get NLConstraint's upper bound
  double GetUpper(const NLConstraint& nlc) const {
    return nlc.GetMainCon().ub();
  }

  /// Get the expression term of an \a NLLogical.
  ExprType GetExpression(const NLLogical& nll) {
    assert( nll.GetResultVar()>=0 );
    return GetInitExpression(nll.GetResultVar());
  }

  /// Get the expression term of an \a NLReification.
  template <int sense>
  ExprType GetExpression(const NLReification<sense>& nll) {
    assert( nll.GetResultVar()>=0 );
    return GetPureInitExpression(nll.GetResultVar());
  }

  /// Get the variable of an \a NLReification.
  template <int sense>
  int GetVariable(const NLReification<sense>& nll) {
    assert( nll.GetResultVar()>=0 );
    return nll.GetResultVar();
  }

  /// GetLinSize(le)
  int GetLinSize(const LinExpression& le) const
  { return le.GetFlatConstraint().GetAffineExpr().size(); }
  /// GetLinCoef(le, i)
  double GetLinCoef(const LinExpression& le, int i) const
  { return le.GetFlatConstraint().GetAffineExpr().coef(i); }
  /// GetLinTerm(le, i)
  Expr GetLinTerm(const LinExpression& le, int i)
  { return GetInitExpression(le.GetFlatConstraint().GetAffineExpr().var(i)); }
  /// GetConstTerm(le)
  double GetConstTerm(const LinExpression& le) const
  { return le.GetFlatConstraint().GetAffineExpr().constant_term(); }

  /// GetLinSize(qe)
  int GetLinSize(const QuadExpression& qe) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetLinTerms().size(); }
  /// GetLinCoef(qe, i)
  double GetLinCoef(const QuadExpression& qe, int i) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetLinTerms().coef(i); }
  /// GetLinTerm(qe, i)
  Expr GetLinTerm(const QuadExpression& qe, int i)
  { return GetInitExpression(qe.GetFlatConstraint().GetQuadExpr().GetBody().GetLinTerms().var(i)); }

  /// GetQuadSize(qe)
  int GetQuadSize(const QuadExpression& qe) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().size(); }
  /// GetQuadCoef(qe, i)
  double GetQuadCoef(const QuadExpression& qe, int i) const
  { return qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().coef(i); }
  /// GetQuadTerm1(qe, i)
  Expr GetQuadTerm1(const QuadExpression& qe, int i)
  { return GetInitExpression(qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().var1(i)); }
  /// GetQuadTerm2(qe, i)
  Expr GetQuadTerm2(const QuadExpression& qe, int i)
  { return GetInitExpression(qe.GetFlatConstraint().GetQuadExpr().GetBody().GetQPTerms().var2(i)); }

  /// GetConstTerm(qe)
  double GetConstTerm(const QuadExpression& qe) const
  { return qe.GetFlatConstraint().GetQuadExpr().constant_term(); }

  /// Get argument expression [\a i]
  template <class FlatExpression>
  Expr GetArgExpression(const FlatExpression& fe, int i)
  { return GetInitExpression(fe.GetFlatConstraint().GetArguments().at(i)); }

  /// Get expression parameter [\a i]
  template <class FlatExpression>
  double GetParameter(const FlatExpression& fe, int i)
  { return fe.GetFlatConstraint().GetParameters().at(i); }

  ////////////////////// INTERNAL ////////////////////////

private:
  /// Get Expr for a variable, if proper/explicit,
  /// or for the InitExpression().
  Expr GetInitExpression(int i_expr) {
    assert(i_expr < is_expr_stored_.size());
    assert(i_expr < expr_stored_.size());
    if (!is_expr_stored_[i_expr]) {
      is_expr_stored_[i_expr] = true;
      if (IsVarProper(i_expr)) {
        expr_stored_[i_expr] = MPD( GetVarExpression(i_expr) );
      } else {
        get_init_expr_(i_expr, &expr_stored_[i_expr]);
      }
    }
    return expr_stored_[i_expr];  // ...............
  }

  /// Get Expr for the InitExpression().
  /// Solver expression for the given variable's init expression.
  /// This is called for NLConstraint
  /// and for NLReification. For them, during result explicification,
  /// the result is marked as variable,
  /// but we still need the expression.
  /// We have to trust that the init expression is accepted.
  /// @note GetInitExpression(\a i_expr) might still return
  ///   the Expr for the result variable \a i_expr.
  Expr GetPureInitExpression(int i_expr) {
    assert(i_expr < is_expr_stored_.size());
    assert(i_expr < expr_stored_.size());
    assert(i_expr < is_init_expr_retrieved_.size());
    assert(!is_init_expr_retrieved_[i_expr]);             // not twice explicified
    if (IsVarProper(i_expr)) {
      is_init_expr_retrieved_[i_expr] = true;             // is being explicified
      Expr result;
      get_init_expr_(i_expr, &result);
      return result;
    }
    return GetInitExpression(i_expr);                     // standard case
  }


public:
  /// Pass vector of var proper flags
  void PassVarProperFlags(std::vector<bool> isvp) {
    is_var_proper_ = std::move(isvp);
    is_expr_stored_.resize(is_var_proper_.size());    // allocate
    expr_stored_.resize(is_var_proper_.size());
    is_init_expr_retrieved_.resize(is_var_proper_.size());
  }

  /// Is var proper?
  bool IsVarProper(int i) const {
    assert(i>=0 && i<(int)is_var_proper_.size());
    return is_var_proper_[i];
  }

  /// Init expr getter type
  using InitExprGetterType = std::function<void(int i_expr, void* pexpr)>;

  /// Provide init expr getter
  void PassInitExprGetter(InitExprGetterType gt)
  { get_init_expr_ = gt; }


private:
  std::vector<bool> is_var_proper_;
  std::vector<bool> is_expr_stored_;    // actual Expr's of the result var or the init expressions
  std::vector<ExprType> expr_stored_;
  std::vector<bool> is_init_expr_retrieved_;
  InitExprGetterType get_init_expr_;
};

}  // namespace mp

#endif // MODEL_API_BASE_H
