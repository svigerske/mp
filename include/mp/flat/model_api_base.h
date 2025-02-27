/*
 Basic flat model API definitions.

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
#ifndef FLAT_MODEL_API_BASE_H_
#define FLAT_MODEL_API_BASE_H_

#include <string>

#include "mp/arrayref.h"
#include "mp/common.h"
#include "mp/flat/constr_base.h"
#include "mp/flat/obj_std.h"
#include "mp/flat/model_info.h"

namespace mp {

/// Define an array of variables
class VarArrayDef {
  ArrayRef<double> lbs_;
  ArrayRef<double> ubs_;
  ArrayRef<var::Type> types_;
  ArrayRef<const char*> names_;
public:
  VarArrayDef() = default;
  template <class BndVec, class TypeVec,
            class NameVec=ArrayRef<const char*> >
  VarArrayDef(BndVec&& lbs, BndVec&& ubs,
              TypeVec&& tys, NameVec&& nms={}) :
    lbs_(std::forward<BndVec>(lbs)),
    ubs_(std::forward<BndVec>(ubs)),
    types_(std::forward<TypeVec>(tys)),
    names_(std::forward<NameVec>(nms)) { }
  VarArrayDef(std::initializer_list<double> lbs,
              std::initializer_list<double> ubs,
              std::initializer_list<var::Type> tys,
              std::initializer_list<const char*> nms = {}) :
    lbs_((lbs)), ubs_((ubs)), types_((tys)), names_(nms) { }
  int size() const { assert(check()); return (int)lbs_.size(); }
  const double* plb() const { return lbs_.data(); }
  const double* pub() const { return ubs_.data(); }
  const var::Type* ptype() const { return types_.data(); }
  void set_lb_ub_types(ArrayRef<double> lbs, ArrayRef<double> ubs,
                       ArrayRef<var::Type> types) {
    lbs_=lbs; ubs_=ubs; types_=types;
    assert(check());
  }
  const char* const* pnames() const { return names_.data(); }
  bool check() const
  { return ubs_.size()==lbs_.size() &&
        types_.size()==lbs_.size(); }
};


/// Level of acceptance of a constraint by a backend
enum class ConstraintAcceptanceLevel {
  NotAccepted=0,
  AcceptedButNotRecommended=1,
  Recommended=2
};

/// Level of acceptance of an expression by a backend
enum class ExpressionAcceptanceLevel {
  NotAccepted=0,
  AcceptedButNotRecommended=1,
  Recommended=2
};


/// ... then for a certain constraint it can be specified
#define ACCEPT_CONSTRAINT(ConstrType, level, con_grp) \
static mp::ConstraintAcceptanceLevel \
    AcceptanceLevel(const ConstrType*) \
{ return mp::ConstraintAcceptanceLevel::level; } \
    static constexpr int \
    GroupNumber(const ConstrType*) { return con_grp; }

/// ... and/or for expressions, e.g., AbsExpression:
#define ACCEPT_EXPRESSION(FlatExprType, level) \
static mp::ExpressionAcceptanceLevel \
AcceptanceLevel(const FlatExprType*) { \
  static_assert( std::is_same_v<FlatExprType, \
      ExprWrapper<typename FlatExprType::FlatConType> >, \
    #FlatExprType \
      " should be an ExprWrapper<> - use standard MP flat expressions" ); \
  return mp::ExpressionAcceptanceLevel::level; \
}        // pc{} checks that FlatExprType = ExprWrapper<...>


/// Constraint group names.
/// @todo Keep consistent with the \a ConstraintGroups enum.
const char* ConGroupName(int cg);

/// Constraint groups
///
/// This is used to access constraint attributes (basis status, duals, ...)
/// Convenient, when the solver accesses constraint attributes in groups.
/// For example, Gurobi 9.5 has linear, quadratic, SOS, and general
///
/// @todo Keep consistent with \a congroup_names
enum ConstraintGroup {
	CG_Default,
	CG_All,
  CG_Algebraic,          // MOSEK 10 seems to have algebraic vs others
  CG_Linear,
  CG_Quadratic,
	CG_Conic,
  CG_General,
  CG_Nonlinear,
  CG_Piecewiselinear,
  CG_SOS,
  CG_SOS1,
  CG_SOS2,
  CG_Logical,
  CG_END_
};


/// ModelAPIs handling custom flat constraints should derive from
class BasicFlatModelAPI {
public:
  /// Placeholder for GetTypeName()
  static const char* GetTypeName()    { return "BasicFlatModelAPI"; }
  /// Placeholder for GetLongName()
  static const char* GetLongName() { return nullptr; }

  /// Placeholder for InitCustomOptions()
  void InitCustomOptions() { }

  /// Pass on a FlatModelInfo object
  void PassFlatModelInfo(const FlatModelInfo* pfmi) {
    pfmi_ = pfmi;
  }

  /// Retrieve FlatModelInfo*
  const FlatModelInfo* GetFlatModelInfo() const { return pfmi_; }

  /// Chance to prepare problem update,
  /// e.g., allocate storage
  void InitProblemModificationPhase(const FlatModelInfo*) {  }
  /// Chance to end problem update
  void FinishProblemModificationPhase() {  }

  ////////////////// Some standard items /////////////////

  /// Placeholder for SetLinearObjective()
  void SetLinearObjective(int , const LinearObjective& ) {
    MP_UNSUPPORTED("FlatModelAPI::SetLinearObjective()");
  }

  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }

  /// Whether accepts NLObjective (relevant in BasicExprModelAPI)
  static int AcceptsNLObj() { return 0; }


  /// Placeholder for SetQuadraticObjective()
  void SetQuadraticObjective(int , const QuadraticObjective& ) {
    MP_UNSUPPORTED("FlatModelAPI::SetQuadraticObjective()");
  }
  /// Placeholder for SetNLObjective()
  void SetNLObjective(int , const NLObjective& ) {
    MP_UNSUPPORTED("FlatModelAPI::SetNLObjective()");
  }

  /// Placeholder for AddConstraint<>()
  template <class Constraint>
  void AddConstraint(const Constraint& ) {
    MP_RAISE(
          std::string("Not handling constraint type '") +
          Constraint::GetTypeName() +
          "'. Provide a handler or a converter method");
  }

  /// Derived backends have to tell C++ to use default handlers if they are needed
  /// when they overload AddConstraint(), due to C++ name hiding
#define USE_BASE_CONSTRAINT_HANDLERS(BaseBackend) \
  using BaseBackend::AddConstraint; \
  using BaseBackend::AcceptanceLevel; \
  using BaseBackend::GroupNumber;

  /// Default constraint group
  static constexpr ConstraintGroup GroupNumber(const BasicConstraint*) {
		return CG_Default;
  }

  /// By default, we say constraint XYZ is not accepted
  static constexpr ConstraintAcceptanceLevel AcceptanceLevel(
      const BasicConstraint*) {
    return ConstraintAcceptanceLevel::NotAccepted;
  }

  /// By default, no expressions (global switch)
  static constexpr ExpressionAcceptanceLevel \
      ExpressionInterfaceAcceptanceLevel()
  { return ExpressionAcceptanceLevel::NotAccepted; }

  /// By default, we say expression XYZ is not accepted
  template <class FlatCon>
  static constexpr ExpressionAcceptanceLevel AcceptanceLevel(
      const ExprWrapper<FlatCon>*) {
    return ExpressionAcceptanceLevel::NotAccepted;
  }

  /// Specifically, ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return false; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin() { return true; }

  /// Specifically, ask if the solver can mix conic quadratic
  /// (entered via dedicated API) and direct quadratic constraints
  static constexpr bool CanMixConicQCAndQC() { return false; }

  /// Ask if the solver can recognize SOCP corner cases
  /// (non-std representations such as xy>=1, see tests)
  /// from quadratic representations
  static constexpr bool CanSOCPCornerCasesFromQC() { return false; }


private:
  const FlatModelInfo* pfmi_ { nullptr };
};


} // namespace mp

#endif // FLAT_MODEL_API_BASE_H_
