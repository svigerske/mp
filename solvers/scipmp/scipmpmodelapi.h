#ifndef SCIPMODELAPI_H
#define SCIPMODELAPI_H

#include "mp/env.h"
#include "scipmpcommon.h"
#include "mp/flat/nl_expr/model_api_base.h"

namespace mp {

class ScipModelAPI :
                     public ScipCommon, public EnvKeeper,
                     public BasicExprModelAPI<ScipModelAPI, SCIP_EXPR*>
{
  using BaseModelAPI = BasicExprModelAPI<ScipModelAPI, SCIP_EXPR*>;

protected:
  void linearHelper(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double lb, const double ub);
  void helpIndicatorLin(const int* pvars, const double* pcoefs, const size_t size, const char* name, const double rhs, const int binary_var, const int binary_value, SCIP_Bool lessthanineq);
  void helpQuad(const char* name, const int* lt_pvars, const double* lt_pcoefs, const size_t lt_size, const int* qt_pvars1, const int* qt_pvars2, const double* qt_pcoefs, const size_t qt_size, const double lb, const double ub);

public:
  /// Construct
  ScipModelAPI(Env& e) : EnvKeeper(e) { }

  /// Class name
  static const char* GetTypeName() { return "ScipModelAPI"; }

  /// If any driver options added from here
  void InitCustomOptions() { }

  /// Called before problem input.
  /// Model info can be used to preallocate memory.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  /// Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
  static int AcceptsQuadObj() { return 0; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  //////////////////////////// EXPRESSION TREES ////////////////////////////
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)
  /// 'NotAccepted' or
  /// 'AcceptedButNotRecommended' would outline each expression
  /// with an auxiliary variable,
  /// but the latter option would enable option acc:_expr,
  /// which when set to 1 switches on the expression trees.
  /// When expressions are off, only flat constraints are used,
  /// e.g., ExpConstraint(a, b), meaning a = exp(b), with a, b variables.
  ///
  /// See also per-expression and per-constraint type switches
  /// (ACCEPT_EXPRESSION and ACCEPT_CONSTRAINT.)
  /// When both are enabled for certain function
  /// (e.g., ExpExpression and ExpConstraint),
  /// option acc:exp allows switching between them.
  ///
  /// Expression trees are submitted to the solver via
  /// the high-level constraints NLConstraint, NLComplementarity,
  /// NLLogical, NLEquivalence, NLImpl, NLRimpl, and NLObjective.
  ACCEPT_EXPRESSION_INTERFACE(AcceptedButNotRecommended);

  /// Whether accepts NLObjective.
  /// No, as of SCIP 9.1.
  static int AcceptsNLObj() { return 0; }

  /// Once expressions are supported, need the following
  /// helper methods.
  ///
  /// GetVarExpression(i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
  SCIP_EXPR* GetVarExpression(int i);

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
  SCIP_EXPR* GetZeroExpression();

  /// For each supported top-level constraint type,
  /// add the ACCEPT_CONSTRAINT macro
  /// and the corresponding AddConstraint function.
  /// Below some typical constraint handlers of a MIP solver.
  /// Further constraint types which could be handled natively by some solvers:
  /// - IndicatorConstraint(Lin/Quad)(LE/EQ/GE)
  /// - Multidirectional indicators Cond(Lin/Quad)Con(LT/LE/EQ/GE/GT), where
  ///   the implication direction (</==/>) depends on the context
  /// - Complementarity
  /// - Logical, counting, piecewise-linear constraints.
  /// See \a constr_std.h.

  /// Below are high-level flat algebraic constraints.

  /// The linear range constraint, if fully supported with basis info etc.
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
  void AddConstraint(const LinConGE& lc);

  /// Ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return true; }

  /// QuadConRange is optional.
  ACCEPT_CONSTRAINT(QuadConRange, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConRange& qc);

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConLE& qc);
  ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConEQ& qc);
  ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
  void AddConstraint(const QuadConGE& qc);

  /// With expression trees, top-level nonlinear algebraic constraints
  /// are submitted to the solver via NLConstraint.
  /// It is still recommended that pure-linear
  /// and pure-quadratic constraints are accepted, then they are
  /// used to submit the corresponding constraint types.
  ///
  /// The implementation can have special treatment
  /// for the case when the linear part has just 1 variable.
  /// This might be so-called variable explicifier var = expr
  /// (or var <= expr, var >= expr.)
  ///
  /// To access information from an NLConstraint,
  /// use the following accessors (don't use methods of NLConstraint itself):
  /// - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
  ///   GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_Nonlinear)
  void AddConstraint(const NLConstraint& nlc);

  /// NL logical constraint: expression = true.
  /// Use GetExpression(nll) to access the expression.
  /// Constraint group: CG_Nonlinear for SCIP,
  /// because using SCIPcreateConsBasicNonlinear().
  ACCEPT_CONSTRAINT(NLLogical, Recommended, CG_Nonlinear)
  void AddConstraint(const NLLogical& nll);

  /// NL equivalence: expression <==> var.
  /// This is an 'expression explicifier'.
  /// Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLEquivalence, Recommended, CG_Nonlinear)
  void AddConstraint(const NLEquivalence& nle);
  /// NL implication: var==1 ==> expression.
  /// This is an 'expression explicifier'.
  /// Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLImpl, Recommended, CG_Nonlinear)
  void AddConstraint(const NLImpl& nle);
  /// NL reverse implication: expression ==> var==1.
  /// This is an 'expression explicifier'.
  /// Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLRimpl, Recommended, CG_Nonlinear)
  void AddConstraint(const NLRimpl& nle);

  /// Moreover, once algebraic expressions are accepted
  /// via NLConstraint, subexpressions might be submitted via
  /// LinExpression and QuadExpression.

  /// LinExpression.
  /// Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(LinExpression, Recommended);
  SCIP_EXPR* AddExpression(const LinExpression& le);

  /// QuadExpression.
  /// Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetQuadSize(le), GetQuadCoef(le, i),
  ///   GetQuadTerm1(le, i), GetQuadTerm2(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(QuadExpression, Recommended);
  SCIP_EXPR* AddExpression(const QuadExpression& le);

  /// Each expression can be accpeted as a proper expression,
  /// or a flat constraint var == expr (with var arguments).
  ///
  /// For each expression,
  /// say ACCEPT_EXPRESSION(Recommended)
  /// and/or ACCEPT_EXPRESSION(AcceptedButNotRecommended).
  /// This can be user-configured via options 'acc:exp' etc.
  ///
  /// Use accessor: GetArgExpression(ee, 0)
  /// - don't ExpExpression's methods.
  ///
  /// Similar for other expression types.
  ACCEPT_EXPRESSION(AbsExpression, Recommended)
  SCIP_EXPR* AddExpression(const AbsExpression& absc);

  /// For each flat constraint type,
  /// say ACCEPT_CONSTRAINT(Recommended)
  /// and/or ACCEPT_CONSTRAINT(AcceptedButNotRecommended).
  /// This can be user-configured via options 'acc:exp' etc.
  ACCEPT_CONSTRAINT(AbsConstraint, Recommended, CG_General)
  void AddConstraint(const AbsConstraint& absc);

  // SCIP 9 has AND/OR as constraints only:
  // ACCEPT_EXPRESSION(AndExpression, AcceptedButNotRecommended)
  // void AddExpression(const AndExpression& cc);
  // ACCEPT_EXPRESSION(OrExpression, Recommended)
  // void AddExpression(const OrExpression& dc);
  ACCEPT_CONSTRAINT(AndConstraint, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const AndConstraint& cc);
  ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  void AddConstraint(const OrConstraint& dc);

  /// Linear indicator constraints can be used as
  /// auxiliary constraints for logical conditions.
  /// If not handled, the compared expressions need
  /// deducible finite bounds for a big-M redefinition.
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

  /// Cones
	ACCEPT_CONSTRAINT(QuadraticConeConstraint, AcceptedButNotRecommended, CG_Conic)
	void AddConstraint(const QuadraticConeConstraint& qc);

  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// Set ``option pl_linearize 0;`` in AMPL if the solver
  /// supports PL natively.
  ACCEPT_CONSTRAINT(SOS1Constraint, Recommended, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, AcceptedButNotRecommended, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);

  /// SCIP nonlinear general constraints and expressions.

  ACCEPT_EXPRESSION(ExpExpression, Recommended)
  SCIP_EXPR* AddExpression(const ExpExpression& );
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended, CG_General)
  void AddConstraint(const ExpConstraint& cc);

  ACCEPT_EXPRESSION(LogExpression, Recommended)
  SCIP_EXPR* AddExpression(const LogExpression& );
  ACCEPT_CONSTRAINT(LogConstraint, Recommended, CG_General)
  void AddConstraint(const LogConstraint& cc);

  /// Use accessor: GetParameter(pe, 0)
  /// - don't use PowExpression's methods.
  ACCEPT_EXPRESSION(PowExpression, Recommended)
  SCIP_EXPR* AddExpression(const PowExpression& );
  ACCEPT_CONSTRAINT(PowConstraint, Recommended, CG_General)
  void AddConstraint(const PowConstraint& cc);

  ACCEPT_EXPRESSION(SinExpression, Recommended)
  SCIP_EXPR* AddExpression(const SinExpression& );
  ACCEPT_CONSTRAINT(SinConstraint, Recommended, CG_General)
  void AddConstraint(const SinConstraint& cc);

  ACCEPT_EXPRESSION(CosExpression, AcceptedButNotRecommended)  //pretty slow in SCIP 8/9
  SCIP_EXPR* AddExpression(const CosExpression& );
  ACCEPT_CONSTRAINT(CosConstraint, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const CosConstraint& cc);

  // TODO Div; PowVarExponent;
  // CondLin... - not really, reader_nl.cpp only handles bool args
};

} // namespace mp

#endif // SCIPMODELAPI_H
