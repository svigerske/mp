#ifndef CONSTR_NL_H
#define CONSTR_NL_H

#include <cmath>
#include <cassert>

#include "mp/error.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Class NLConstraint.
/// Algebraic range constraint with a linear part
/// and an expression term: `lb <= a'x + expr <= ub`.
/// LinConRange is a member to avoid overloading
/// when deriving from an existing constraint type.
class NLConstraint
    : public BasicConstraint, public NumericFunctionalConstraintTraits {
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
      : lcr_(lt, rng, fSort), expr_(expr) { }

  /// Get the main algebraic constraint
  const LinConRange& GetMainCon() const { return lcr_; }

  /// Has expression term?
  bool HasExpr() const { return expr_>=0; }

  /// Expression index.
  /// @note ModelAPI should call self.HasExpression()
  ///   and self.GetExpression()
  ///   to obtain the expression term.
  int ExprIndex() const { assert(HasExpr()); return expr_; }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  /// Compute violation.
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x, bool logical=false) const {
    double bd = lcr_.GetBody().ComputeValue(x);
    if (HasExpr())
      bd += x[ExprIndex()];     // Add expr value. Assume it's precomputed
    if (!logical) {
      if (lcr_.lb() > bd)
        return {lcr_.lb() - bd, lcr_.lb()};
      if (bd > lcr_.ub())
        return {bd - lcr_.ub(), lcr_.ub()};
      return
          {std::max( // negative. Same for strict cmp?
               lcr_.lb() - bd, bd - lcr_.ub()),
           0.0};
    }
    return {double(!lcr_.is_valid(bd)), 1.0};
  }


private:
  LinConRange lcr_;
  int expr_ {-1};
};


/// Export to JSON
inline void WriteJSON(JSONW jw,
                      const NLConstraint& nlc) {
  WriteJSON(jw["main_alg_con"], nlc.GetMainCon());
  if (nlc.HasExpr())
    WriteJSON(jw["expr_index"], nlc.ExprIndex());
}

/// Write RhsCon without name.
template <class Writer, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLConstraint& nlc,
                           const Names& vnam) {
  wrt << "NLExprIndex: " << vnam.at(nlc.ExprIndex()) << " IN: ";
  WriteModelItem(wrt, nlc.GetMainCon(), vnam);
}


/// Algebraic expression explicifier.
/// Syntax sugar for the assignment: var <=/==/>= expr.
/// Can have special meaning in certain solvers.
/// Sense: equality (0), >= (1), <= (-1).
/// Can be implemented as == for all senses (e.g., GRBaddgenconstrNL),
/// the inequalities can be used to preserve convexity.
/// This is a static constraint.
template <int sense>
class NLBaseAssign
    : public BasicConstraint, public NumericFunctionalConstraintTraits {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    if (0==sense) return "NLAssignEQ";
    if (-1==sense) return "NLAssignLE";
    if (1==sense) return "NLAssignGE";
    MP_RAISE("NLBaseAssign: unknown sense");
  }

  /// Construct
  NLBaseAssign(int b) : bvar_(b) { }

  /// Get var
  int GetVar() const { return bvar_; }

  /// Imitate an algebraic con
  int size() const { return 1; }
  /// Imitate
  double coef(int i) const { assert(!i); return -1.0; }
  /// Imitate
  double var(int i) const { assert(!i); return GetVar(); }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  // Compute violation... Should be 0

private:
  int bvar_ {-1};
};


/// Typedef NLAssignEQ
using NLAssignEQ = NLBaseAssign<0>;
/// Typedef NLAssignLE
using NLAssignLE = NLBaseAssign<-1>;
/// Typedef NLAssignGE
using NLAssignGE = NLBaseAssign<1>;

/// Write a Reification
template <int sense>
inline void WriteJSON(JSONW jw,
                      const NLBaseAssign<sense>& reif) {
  jw["var_explicit_assign"] = reif.GetVar();
  jw["sense"] = sense;
}

/// Write RhsCon without name.
template <class Writer, int sense, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLBaseAssign<sense>& nlr,
                           const Names& vnam) {
  wrt << "EXPLICIT ASSIGN var: " << vnam.at(nlr.GetVar());
}


/// NLComplementarity
/// TODO extra class, to enable ACCEPT_CONSTRAINT
using NLComplementarity = ComplementarityConstraint<AffineExpr>;


/// NL logical constraint: expr(resvar) == true
class NLLogical
    : public BasicConstraint {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    return "NLLogical";
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Construct from the result variable
  NLLogical(int rv) : resvar_(rv) { assert(rv>=0); }

  /// Get resvar
  int GetCapturedResultVar() const { return resvar_; }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  /// Compute violation
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x) const {
    return {std::fabs(x[GetCapturedResultVar()] - 1.0), 1.0};
  }

private:
  int resvar_ {-1};
};

/// Write an NLLogical
inline void WriteJSON(JSONW jw,
                      const NLLogical& nll) {
  jw["resvar"] = nll.GetCapturedResultVar();
}

/// Write RhsCon without name.
template <class Writer, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLLogical& nllc,
                           const Names& vnam) {
  wrt << "NLLogicalExprIndex: "
      << vnam.at(nllc.GetCapturedResultVar());
}


/// NL equivalence constraint: expr(var1) <==> expr(var2).
/// We have to do this new kind because flat constraints
/// represent equivalence algebraically.
class NLEquivalence
    : public BasicConstraint {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    return "NLEquivalence";
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Construct from the result variable
  NLEquivalence(int v1, int v2) : v1_(v1), v2_(v2)
  { assert(v1>=0 && v2>=0); }

  /// Get var 1
  int GetVar1() const { return v1_; }

  /// Get var 2
  int GetVar2() const { return v2_; }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  /// Compute violation
  template <class VarInfo>
  Violation
  ComputeViolation(const VarInfo& x) const {
    return {std::fabs(x[GetVar1()] - x[GetVar2()]), 1.0};
  }

private:
  int v1_ {-1}, v2_ {-1};
};

/// Write an NLLogical
inline void WriteJSON(JSONW jw,
                      const NLEquivalence& nll) {
  jw["var1"] = nll.GetVar1();
  jw["var2"] = nll.GetVar2();
}

/// Write RhsCon without name.
template <class Writer, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLEquivalence& nllc,
                           const Names& vnam) {
  wrt << "NLEquivalenceVar1: " << vnam.at(nllc.GetVar1());
  wrt << "NLEquivalenceVar2: " << vnam.at(nllc.GetVar2());
}


/// Logical expression explicifier.
/// Syntax sugar for reification: b==1 <==> expr(b)==1.
/// Sense: equivalence (0), impl(-1), rimpl (1).
/// This is a static constraint.
template <int sense>
class NLBaseReif
    : public BasicConstraint {
public:
  /// Constraint type name
  static const char* GetTypeName() {
    if (0==sense) return "NLReifEquiv";
    if (-1==sense) return "NLReifImpl";
    if (1==sense) return "NLReifRimpl";
    MP_RAISE("NLReif: unknown sense");
  }

  /// Is logical?
  static bool IsLogical() { return true; }

  /// Construct
  NLBaseReif(int b) : bvar_(b) { }

  /// Get bvar
  int GetBVar() const { return bvar_; }

  /// Throw - should not be used
  VarArray1 GetArguments() const { MP_RAISE("No marking for NL items"); }

  // Compute violation... Should be 0

private:
  int bvar_ {-1};
};


/// Typedef NLReifEquiv
using NLReifEquiv = NLBaseReif<0>;
/// Typedef NLReifImpl
using NLReifImpl = NLBaseReif<-1>;
/// Typedef NLReifRImpl
using NLReifRimpl = NLBaseReif<1>;

/// Write a Reification
template <int sense>
inline void WriteJSON(JSONW jw,
                      const NLBaseReif<sense>& reif) {
  jw["var_explicit_reif"] = reif.GetBVar();
  jw["sense"] = sense;
}

/// Write RhsCon without name.
template <class Writer, int sense, class Names>
inline void WriteModelItem(Writer& wrt,
                           const NLBaseReif<sense>& nlr,
                           const Names& vnam) {
  wrt << "EXPLICIT REIF var: " << vnam.at(nlr.GetBVar());
}

}  // namespace mp

#endif // CONSTR_NL_H
