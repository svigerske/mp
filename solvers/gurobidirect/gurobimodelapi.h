#ifndef GUROBIMODELAPI_H
#define GUROBIMODELAPI_H

#include <memory>

#include "mp/env.h"
#include "gurobicommon.h"
#include "mp/flat/std_constr.h"

namespace mp {


class GurobiModelAPI : public GurobiCommon, public EnvKeeper {
  using BaseModelAPI = BasicFlatModelAPI;

public:
  /// Model API name
  static const char* GetBackendName();
  /// Unused
  static const char* GetBackendLongName() { return nullptr; }

  /// Construct
  GurobiModelAPI(Env& e) : EnvKeeper(e) { }

  /// This is called before model is pushed to the Backend
  void InitProblemModificationPhase();
  /// Chance to call GRBupdatemodel()
  void FinishProblemModificationPhase();

  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////// MODELING ACCESSORS /////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective( int iobj, const QuadraticObjective& qo );
  /// First objective's sense
  void NoteGurobiMainObjSense(obj::Type s);
  obj::Type GetGurobiMainObjSense() const;
  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  /// TODO Linear constraint attributes (lazy/user cut, etc)
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);
  ACCEPT_CONSTRAINT(QuadraticConstraint, Recommended, CG_Quadratic)
  void AddConstraint(const QuadraticConstraint& qc);
  ACCEPT_CONSTRAINT(MaximumConstraint, acc_max(), CG_General)
  void AddConstraint(const MaximumConstraint& mc);
  ACCEPT_CONSTRAINT(MinimumConstraint, acc_min(), CG_General)
  void AddConstraint(const MinimumConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, acc_abs(), CG_General)
  void AddConstraint(const AbsConstraint& absc);
  ACCEPT_CONSTRAINT(ConjunctionConstraint, acc_and(), CG_General)
  void AddConstraint(const ConjunctionConstraint& cc);
  ACCEPT_CONSTRAINT(DisjunctionConstraint, acc_or(), CG_General)
  void AddConstraint(const DisjunctionConstraint& mc);
  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, acc_ind_le(), CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, acc_ind_eq(), CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);

  /// General
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_General)
  void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, Recommended, CG_General)
  void AddConstraint(const ExpAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_General)
  void AddConstraint(const LogConstraint& cc);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended, CG_General)
  void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
  void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended, CG_General) // y = cos(x)
  void AddConstraint(const CosConstraint& cc);  // GRBaddgenconstrCos(x, y);
  ACCEPT_CONSTRAINT(TanConstraint, Recommended, CG_General)
  void AddConstraint(const TanConstraint& cc);
  ACCEPT_CONSTRAINT(PLConstraint, Recommended, CG_General)
  void AddConstraint(const PLConstraint& cc);

  void InitOptions();

private:
  /// The sense of the main objective
  obj::Type main_obj_sense_;

  /// These options are stored in the class as variables
  /// for direct access
  struct Options {
    int acc_min_=2, acc_max_=2, acc_abs_=2, acc_and_=2, acc_or_=2,
      acc_ind_le_=2, acc_ind_eq_=2;
  } storedOptions_;


protected:  //////////// Option accessors ////////////////
  int acc_abs() const { return storedOptions_.acc_abs_; }
  int acc_min() const { return storedOptions_.acc_min_; }
  int acc_max() const { return storedOptions_.acc_max_; }
  int acc_and() const { return storedOptions_.acc_and_; }
  int acc_or() const { return storedOptions_.acc_or_; }
  int acc_ind_le() const { return storedOptions_.acc_ind_le_; }
  int acc_ind_eq() const { return storedOptions_.acc_ind_eq_; }
};

} // namespace mp

#endif // GUROBIMODELAPI_H
