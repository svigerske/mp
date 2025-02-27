#ifndef MOSEKMODELAPI_H
#define MOSEKMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "mosekcommon.h"
#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

class MosekModelAPI :
    public MosekCommon, public EnvKeeper,
    public BasicFlatModelAPI
{
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Construct
  MosekModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "MosekModelAPI"; }

  /// Called before problem modification.
  /// @param fmi: current problem information
  void InitProblemModificationPhase(const FlatModelInfo* fmi);
  /// After
  void FinishProblemModificationPhase();

  /// TODO Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 1; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// TODO For each supported constraint type, add the ACCEPT_CONSTRAINT macro
  /// and the relative AddConstraint function.
  /// Below some typical constraint handlers of a MIP solver.
  /// Further constraint types which could be handled natively by some solvers:
  /// - IndicatorConstraint(Lin/Quad)(LE/EQ/GE)
  /// - Multidirectional indicators Cond(Lin/Quad)Con(LT/LE/EQ/GE/GT), where
  ///   the implication direction (</==/>) depends in the expression context
  /// - Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constr_std.h and other drivers.


  /// The linear range constraint, if fully supported with basis info etc.
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Algebraic)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Algebraic)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Algebraic)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Algebraic)
  void AddConstraint(const LinConGE& lc);

  /// QuadConRange is optional.
  ACCEPT_CONSTRAINT(QuadConRange, Recommended, CG_Algebraic)
  void AddConstraint(const QuadConRange& qc);

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Algebraic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Algebraic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Algebraic)
  void AddConstraint(const QuadConGE& qc);

  /// Specifically, ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return false; }

  /// If cvt:prod=7 (and not 5) default.
  /// Recommendation to return the opposite value as
  /// AcceptsNonconvexQC().
  static constexpr bool WantLogicalizedProd2Bin()
  { return !AcceptsNonconvexQC(); }

  /// Specifically, ask if the solver can mix conic quadratic
  /// (entered via dedicated API) and direct quadratic constraints
  static constexpr bool CanMixConicQCAndQC() { return false; }

  /// Cones
  ACCEPT_CONSTRAINT(QuadraticConeConstraint, Recommended, CG_Conic)
  void AddConstraint(const QuadraticConeConstraint& qc);
  ACCEPT_CONSTRAINT(RotatedQuadraticConeConstraint, Recommended, CG_Conic)
  void AddConstraint(const RotatedQuadraticConeConstraint& qc);
  ACCEPT_CONSTRAINT(ExponentialConeConstraint, Recommended, CG_Conic)
  void AddConstraint(const ExponentialConeConstraint& ec);

  /// Linear indicator constraints can be used as
  /// auxiliary constraints for logical conditions.
  /// Mosek 10 adds DJC which we use for indicators.
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// Set ``option pl_linearize 0;`` in AMPL if the solver
  /// supports PL natively.
  /// MOSEK 10 has no SOS
  ACCEPT_CONSTRAINT(SOS1Constraint, NotAccepted, CG_General)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, NotAccepted, CG_General)
  void AddConstraint(const SOS2Constraint& cc);


protected:
  /// internal, add linear constraint
  void AddLinearConstraint(
      MSKtask_t lp, size_t size, MSKboundkey_enum key,
      double lb, double ub,
      const int* vindex, const double* values, const char* name);


private:
  int n_alg_cons_ = 0;
};

} // namespace mp

#endif // MOSEKMODELAPI_H
