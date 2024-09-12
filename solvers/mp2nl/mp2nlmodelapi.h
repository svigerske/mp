#ifndef MP2NLMODELAPI_H
#define MP2NLMODELAPI_H

#include "mp/env.h"
#include "mp2nlcommon.h"
#include "mp/flat/nl_expr/model_api_base.h"

#include "mp/nl-feeder.h"


namespace mp {

/// Expression ID for MP2NLModelAPI.
/// Can be
/// 1. "Empty expression" (0.0),
/// 2. variable expression (represent a variable),
/// 3. normal expression node.
class MP2NL_Expr {
public:
  /// Construct
  MP2NL_Expr(int e=0) : id_(e) { }

	/// Get the expr ID
	int GetID() const { return id_; }

private:
	int id_ {0};
};


/// MP2NLModelAPI.
/// MP2NL is to be used as a meta-driver,
/// performing reformulations for the final NL solver.
///
/// MP2NLModelAPI translates the reformulated model for the NL Writer.
class MP2NLModelAPI
		: public MP2NLCommon, public EnvKeeper,
      public BasicExprModelAPI<MP2NLModelAPI, MP2NL_Expr>,
      public NLFeeder<MP2NLModelAPI, MP2NL_Expr>
{
public:
  using BaseModelAPI = BasicExprModelAPI<MP2NLModelAPI, MP2NL_Expr>;
	using Expr = MP2NL_Expr;
  using BaseNLFeeder = NLFeeder<MP2NLModelAPI, MP2NL_Expr>;

  /// Construct
  MP2NLModelAPI(Env& e) : EnvKeeper(e) {
    CreateInterfaces();
  }

  /// Class name
	static const char* GetTypeName() { return "MP2NLModelAPI"; }

  /// If any driver options added from the ModelAPI
  /// (usually only in the Backend)
  void InitCustomOptions() { }

  /// Called before problem input.
  /// Model info can be used to preallocate memory.
  void InitProblemModificationPhase(const FlatModelInfo*);
  /// After
  void FinishProblemModificationPhase();

  /// Implement the following functions using the solver's API
  void AddVariables(const VarArrayDef& vad);
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  /// Whether accepting quadratic objectives:
  /// 0 - no, 1 - convex, 2 - nonconvex
	static int AcceptsQuadObj() { return 2; }
  void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);
  /// Whether accepts NLObjective.
  static int AcceptsNLObj() { return 1; }
  void SetNLObjective(int iobj, const NLObjective& nlo);


  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  /// Handle flat constraints: inherit basic API
  USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

  //////////////////////////// EXPRESSION TREES ////////////////////////////
  /// Hanlde expression trees: inherit basic API
  USE_BASE_EXPRESSION_HANDLERS(BaseModelAPI)


  ///////////// !!! First go with linear models. //////////////////


  /// ACCEPT_EXPRESSION_INTERFACE():
  /// 'NotAccepted' or
  /// 'AcceptedButNotRecommended' would resort to flat constraints uniformly
  /// which outline each expression with an auxiliary variable:
  /// e.g., ExpConstraint(a, b), meaning a = exp(b), with a, b variables.
  /// But the latter option would enable the solver option acc:_expr,
  /// which, when set to 1, switches on the expression trees.
  ///
  /// See also per-expression and per-constraint type switches
  /// (ACCEPT_EXPRESSION and ACCEPT_CONSTRAINT.)
  /// When both are enabled for certain function
  /// (e.g., CosExpression and CosConstraint),
  /// solver option acc:cos allows switching between them.
  ///
  /// Expression trees are submitted to the solver via
  /// the high-level constraints
  /// NLConstraint, NLAssignLE, NLAssignEQ, NLAssignGE,
  /// NLComplementarity,
  /// NLLogical, NLEquivalence, NLImpl, NLRimpl, and NLObjective.
  ACCEPT_EXPRESSION_INTERFACE(AcceptedButNotRecommended);

  /// Once expressions are supported, need the following
  /// helper methods.
  ///
  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
	Expr GetVarExpression(int i);

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
	Expr GetZeroExpression();

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
	/// NL format supports this.
  ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Algebraic)
  void AddConstraint(const LinConRange& lc);

  /// LinCon(LE/EQ/GE) should have 'Recommended' for all non-expression backends
  /// and have an implementation,
  /// or a conversion rule is needed in a derived FlatConverter
  /// Even for expression backends, they can be implemented for efficiency.
  ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Algebraic)
  void AddConstraint(const LinConLE& lc);
  ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Algebraic)
  void AddConstraint(const LinConEQ& lc);
  ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Algebraic)
  void AddConstraint(const LinConGE& lc);

  /// Ask if the solver accepts non-convex quadratic constraints
  static constexpr bool AcceptsNonconvexQC() { return true; }

  /// QuadConRange is optional.
	/// No pure variable-quadratics for NL format.
	ACCEPT_CONSTRAINT(QuadConRange, NotAccepted, CG_Quadratic)

  /// If using quadratics,
  /// QuadCon(LE/EQ/GE) should have 'Recommended'
  /// and have an implementation.
  /// Even for expression backends, they can be implemented for efficiency,
  /// if the solver supports pure variable-quadratics.
  ACCEPT_CONSTRAINT(QuadConLE, NotAccepted, CG_Quadratic)
	ACCEPT_CONSTRAINT(QuadConEQ, NotAccepted, CG_Quadratic)
	ACCEPT_CONSTRAINT(QuadConGE, NotAccepted, CG_Quadratic)


  /// Cones
  ACCEPT_CONSTRAINT(QuadraticConeConstraint, NotAccepted, CG_Conic)


  /// If NLConstraint is accepted, with expression trees,
  /// then top-level nonlinear algebraic constraints
  /// are submitted to the solver via NLConstraint.
  ///
  /// @note If NLConstraint is not accepted, then LinCon(LE/GE/EQ/[Range])
  ///   should be. In such case, top-level algebraic expressions
  ///   are explicified via NLAssign(LE/GE/EQ),
  ///   for example, for Gurobi 12.
  ///
  /// @note It is still recommended that pure-linear
  ///   and pure-quadratic constraints are accepted, then they are
  ///   used to submit the corresponding constraint types
  ///   without expressions.
  ///
  /// @note To access information from an NLConstraint,
  ///   use the following accessors (don't use methods of NLConstraint itself):
  ///   - GetLinSize(nlc), GetLinCoef(nlc, i), GetLinVar(nlc, i),
  ///     GetExpression(nlc), GetLower(nlc), GetUpper(nlc).
  ACCEPT_CONSTRAINT(NLConstraint, Recommended, CG_Algebraic)
  void AddConstraint(const NLConstraint& nlc);

  /// NLAssignEQ: algebraic expression expicifier.
  /// Meaning: var == expr.
  ///
  /// @note Should be 'Recommended' for expression solvers.
  /// @note Example: GRBaddgenconstrNL in Gurobi 12.
  ///   In other API types can be implemented using
  ///   NL algebraic constraint.
  ///
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignEQ, Recommended, CG_Algebraic)
  void AddConstraint(const NLAssignEQ& nle);
  /// NLAssignLE: algebraic expression expicifier in positive context.
  /// Meaning: var <= expr.
  ///
  /// @note Should be 'Recommended' in expression solvers.
  /// @note Can be implemented as NLAssignEQ,
  ///   but this may lose convexity.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignLE, Recommended, CG_Algebraic)
  void AddConstraint(const NLAssignLE& nle);
  /// NLAssignGE: algebraic expression expicifier in negative context.
  /// Meaning: var >= expr.
  ///
  /// @note Should be 'Recommended' in expression solvers.
  /// @note Can be implemented as NLAssignEQ,
  ///   but this may lose convexity.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLAssignGE, Recommended, CG_Algebraic)
  void AddConstraint(const NLAssignGE& nle);

  /// @todo
  ACCEPT_CONSTRAINT(NLComplementarity, NotAccepted, CG_Algebraic)
  void AddConstraint(const NLComplementarity& cc);


  /// NL logical constraint: expression = true.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Use GetExpression(nll) to access the expression.
  ACCEPT_CONSTRAINT(NLLogical, Recommended, CG_Logical)
  void AddConstraint(const NLLogical& nll);

  /// NL equivalence: expression <==> var.
  /// This is an 'expression explicifier'.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLEquivalence, Recommended, CG_Logical)
  void AddConstraint(const NLEquivalence& nle);
  /// NL implication: var==1 ==> expression.
  /// This is an expression explicifier in positive context.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLImpl, Recommended, CG_Logical)
  void AddConstraint(const NLImpl& nle);
  /// NL reverse implication: expression ==> var==1.
  /// This is an expression explicifier in negative context.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLRimpl, Recommended, CG_Logical)
  void AddConstraint(const NLRimpl& nle);


  /// Linear indicator constraints can be used as
  /// auxiliary constraints for logical conditions.
  /// If not handled, the compared expressions need
  /// deducible finite bounds for a big-M redefinition.
  ///
  /// @note Need them in NL format output? Currently seems yes,
  ///   several reformulations produce it.
  ///   As long as the target solver accepts.
  /// @todo But if ever there appears a solver that only accepts
  ///   Cond(Lin/Quad)..., we should produce them instead.
  /// Mosek???
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, NotAccepted, CG_Logical)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, NotAccepted, CG_Logical)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, NotAccepted, CG_Logical)
  void AddConstraint(const IndicatorConstraintLinGE& mc);


  /// SOS constraints can be used by AMPL for redefinition of
  /// piecewise-linear expressions.
  /// @note Set ``option pl_linearize 0;`` in AMPL if the solver
  ///   supports PL natively.
  ACCEPT_CONSTRAINT(SOS1Constraint, NotAccepted, CG_SOS)
  void AddConstraint(const SOS1Constraint& cc);
  ACCEPT_CONSTRAINT(SOS2Constraint, NotAccepted, CG_SOS)
  void AddConstraint(const SOS2Constraint& cc);



  /// @brief LinExpression.
  /// Once algebraic expressions are accepted
  /// via NLConstraint, subexpressions might be submitted via
  /// LinExpression and QuadExpression.
  ///
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(LinExpression, Recommended);
  Expr AddExpression(const LinExpression& le);

  /// QuadExpression.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetQuadSize(le), GetQuadCoef(le, i),
  ///   GetQuadTerm1(le, i), GetQuadTerm2(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(QuadExpression, Recommended);
  Expr AddExpression(const QuadExpression& le);

  /// Each expression can be accepted as a proper expression,
  /// or as a flat functional constraint var <=/==/>= expr
  /// (in this case, with variables as arguments).
  /// The uequality/inqeuality type of the flat constraint is
  /// determied by GetContext().
  ///
  /// @note For each expression,
  /// say ACCEPT_EXPRESSION(Recommended)
  /// and/or ACCEPT_EXPRESSION(AcceptedButNotRecommended).
  /// This can be user-configured via solver options 'acc:exp' etc.
  ///
  /// @note Use accessor: GetArgExpression(ee, 0)
  /// - don't ExpExpression's methods.
  ///
  /// Similar for other expression types.
  ACCEPT_EXPRESSION(AbsExpression, Recommended)
  Expr AddExpression(const AbsExpression& absc);

  /// For each flat constraint type,
  /// say ACCEPT_CONSTRAINT(Recommended)
  /// and/or ACCEPT_CONSTRAINT(AcceptedButNotRecommended).
  /// This can be user-configured via solver options 'acc:exp' etc.
  ACCEPT_CONSTRAINT(AbsConstraint, NotAccepted, CG_General)

  ACCEPT_EXPRESSION(AndExpression, Recommended)
  MP2NL_Expr AddExpression(const AndExpression& cc);
  ACCEPT_EXPRESSION(OrExpression, Recommended)
  MP2NL_Expr AddExpression(const OrExpression& dc);
  // ACCEPT_CONSTRAINT(AndConstraint, AcceptedButNotRecommended, CG_General)
  // void AddConstraint(const AndConstraint& cc);
  // ACCEPT_CONSTRAINT(OrConstraint, Recommended, CG_General)
  // void AddConstraint(const OrConstraint& dc);


  ACCEPT_EXPRESSION(ExpExpression, Recommended)
  Expr AddExpression(const ExpExpression& );
  ACCEPT_EXPRESSION(LogExpression, Recommended)
  Expr AddExpression(const LogExpression& );
  /// @note Use accessor: GetParameter(pe, 0)
  ///   - don't use PowExpression's methods.
  ACCEPT_EXPRESSION(PowExpression, Recommended)
  Expr AddExpression(const PowExpression& );
  ACCEPT_EXPRESSION(SinExpression, Recommended)
  Expr AddExpression(const SinExpression& );
  ACCEPT_EXPRESSION(CosExpression, Recommended)
  Expr AddExpression(const CosExpression& );

  // TODO Div; PowVarVar;
  // CondLin... - not really, reader_nl.cpp only handles bool args


public:
  ///////////////////////////////////////////////////////////////////
  ///////////////// Implement NLFeeder interface ////////////////////
  ///////////////////////////////////////////////////////////////////

  ///////////////////// 1. NL HEADER AND OPTIONS /////////////////
  /** Provide NLHeader.
     *
     *	This method is called first.
     *
     *  NLHeader summarizes the model and provides some
     *  technical parameters,
     *  such as text/binary NL format. */
  NLHeader Header();

  /// NL comments?
  bool WantNLComments() const { return false; }

  /// The maximum number of significant digits written.
  /// The default value requests full precision, which
  /// might be the shortest representation that, when
  /// converted to binary and properly rounded, will
  /// give exactly the binary value stored in the computer.
  int OutputPrecision() const { return 0; }

  /// Write bounds first?
  /// The default is yes in AMPL, controlled by
  /// (the value of option nl_permute) & 32
  /// (the bit is 0 for yes).
  /// Changing this option is deprecated, see
  /// https://netlib.org/ampl/changes.
  bool WantBoundsFirst() const { return true; }

  /// Want Jacobian column sizes?
  /// Required by some nonlinear solvers.
  /// Options: 0 - none, 1 - cumulative,
  /// 2 - non-cumulative.
  /// This option controls how ColSizeWriter
  /// writes the provided sizes (which should be
  /// non-cumulative).
  int WantColumnSizes() const { return 1; }


  ///////////////////// 2. OBJECTIVES /////////////////////
  /** Description for objective function \a i
   *    (\a i in 0..num_objs-1).
   *  With WantNLComments()==true, this is
     *  written to text-format NL as a comment. */
  const char* ObjDescription(int i) { return ""; }

  /** Provide type of objective \a i.
     *  - 0 - minimization;
     *  - 1 - maximization. */
  int ObjType(int i) {
    switch (obj_info_[i].GetStaticTypeID()) {
    case StaticItemTypeID::ID_LinearObjective:
    case StaticItemTypeID::ID_NLObjective:
      return obj::MIN==((LinearObjective*)(obj_info_[i].GetPItem()))->obj_sense()
          ? 0 : 1;
    default:
      MP_RAISE("Unknown objective type");
    }
  }

  /** Feed gradient for objective \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern).
   *
   *  Implementation skeleton:
   *      if (obj_grad[i].size()) {
   *        auto svw = svwf.MakeVectorWriter(obj_grad[i].size());
   *        for (size_t j=0; j<obj_grad.size(); ++j)
   *          svw.Write(obj_grad[j].var_index, obj_grad[j].coef);
   *      }
   */
  template <class ObjGradWriterFactory>
  void FeedObjGradient(int i, ObjGradWriterFactory& svwf) {
    switch (obj_info_[i].GetStaticTypeID()) {
    case StaticItemTypeID::ID_LinearObjective:
    case StaticItemTypeID::ID_NLObjective: {
      const auto& obj = *((LinearObjective*)(obj_info_[i].GetPItem()));
      if (obj.num_terms()) {
        auto svw = svwf.MakeVectorWriter(obj.num_terms());
        for (int j=0; j<obj.num_terms(); ++j)
          svw.Write(GetNewVarIndex( obj.GetLinTerms().var(j) ),      // new ordering
                    obj.GetLinTerms().coef(j));
      }
    } break;
    default:
      MP_RAISE("Unknown objective type");
    }
  }

  /** Feed nonlinear expression of objective \a i.
   *
   *  The default implementation below feeds constant 0
   *  (linear models.)
   *
   *  Implementation example:
   *      ew.EPut(obj_root_expr[i]);
   *
   *  Details of ObjExprWriter: see NLWriter2. */
  template <class ObjExprWriter>
  void FeedObjExpression(int , ObjExprWriter& ew)
  { ew.NPut(0.0); }


  ///////////////////// 3. DEFINED VARIABLES /////////////////////
  /** Defined variables.
     *
     *  Classical NL writes first the defined variables
     *  which are used in several places (constraints and/or
     *  objectives). Defined variables used in a single place
     *  (1 constraint, or 1 objective), are written
     *  just before the expression tree of their usage.
     *
     *  For most solvers, this requirement can be ignored
     *  and this method can return all defined variables
     *  in the first group (for \a i=0).
     *
     *	The method is guaranteed to be called in the following order:
     *		1. For \a i=0;
     *		2. For \a i>0, increasing, before constraint \a (i-1)'s expression;
     *		3. For \a i<0, decreasing, before objective \a (-i-1)'s expression.
     *
     *  @param i:
     *		- For \a i=0, feed a sequence of defined variables
     *			used in several constraints and/or objectives.
     *		- For \a i>0, feed the defined variables used solely
     *			in constraint \a i-1.
     *		- For \a i<0, feed the defined variables used solely
     *			in objective \a -i-1.
     *
   *  Implementation skeleton:
   *      // dvar_index in num_vars..num_vars+num_defvars-1.
   *      for (int dvar_index: dvar_indexes[i]) {
   *        auto dv = dvw.StartDefVar(dvar_index, lin_nnz, name_or_comment);
   *        /////////// Write the linear part:
   *        auto linw = dv.GetLinExprWriter();
   *        for (int i=0; i<lin_nnz; ++i)
   *          linw.Write(linexp_var[i], linexp_coef[i]);
   *        /////////// Write the expression tree:
   *        auto ew = dv.GetExprWriter();
   *        ew.EPut(root_expr);
   *      }
     */
  template <class DefVarWriterFactory>
  void FeedDefinedVariables(int i, DefVarWriterFactory& ) { }


  ///////////////////// 4. VARIABLE BOUNDS /////////////////////
  /** Bounds for variables (except defined variables).
   *  Use +-inf for missing lower and/or upper bounds.
     *  Note that variable type is given by variable ordering,
   *  see NLHeader.
   *
   *  Implementation skeleton:
   *      for (int i = 0; i < hdr.num_vars; i++)
   *        vbw.WriteLbUb(lb[i], ub[i]);
   */
  template <class VarBoundsWriter>
  void FeedVarBounds(VarBoundsWriter& vbw) {
    for (size_t i = 0; i < var_lbs_.size(); i++) {
      auto i_old = GetOldVarIndex(i);
      vbw.WriteLbUb(var_lbs_[i_old], var_ubs_[i_old]);
    }
  }


  ///////////////// 5. CONSTRAINT BOUNDS & COMPLEMENTARITY ///////
  /// \rst
  /// Algebraic constraint bounds (for a single constraint):
  /// either range (lb, ub),
  /// or complementarity info (k, cvar), when k>0.
  ///
  /// For a complementarity constraint to hold, if cvar is at
  ///	its lower bound, then body >= 0; if cvar is at its upper
  /// bound, then body <= 0;
  ///	and if cvar is strictly between its bounds, then body = 0.
  /// The integer k in a complementarity constraint line indicates
  /// which bounds on cvar are finite: 1 and 3 imply a finite
  /// lower bound; 2 and 3 imply a finite upper bound; 0 (which
  ///	should not occur) would imply no finite bounds, i.e.,
  /// body = 0 must always hold.
  ///
  /// Example:
  ///
  /// .. code-block:: ampl
  ///
  ///    ampl: var x; var y; var z;
  ///	   ampl: s.t. Compl1: x+y >= 3 complements x-z <= 15;
  ///	   ampl: s.t. Compl2: -2 <= 2*y+3*z <= 13 complements 6*z-2*x;
  ///	   ampl: expand;
  ///	   subject to Compl1:
  ///					3 <= x + y
  ///			 complements
  ///					x - z <= 15;
  ///
  ///	   subject to Compl2:
  ///					-2 <= 2*y + 3*z <= 13
  ///			 complements
  ///					-2*x + 6*z;
  ///
  ///	   ampl: solexpand;
  ///	   Nonsquare complementarity system:
  ///					4 complementarities including 2 equations
  ///					5 variables
  ///	   subject to Compl1.L:
  ///					x + y + Compl1$cvar = 0;
  ///
  ///	   subject to Compl1.R:
  ///					-15 + x - z <= 0
  ///			 complements
  ///					Compl1$cvar <= -3;
  ///
  ///	   subject to Compl2.L:
  ///					2*y + 3*z - Compl2$cvar = 0;
  ///
  ///	   subject to Compl2.R:
  ///					-2*x + 6*z
  ///			 complements
  ///					-2 <= Compl2$cvar <= 13;
  ///
  /// \endrst
  using BaseNLFeeder::AlgConRange;

  /** Bounds/complementarity for all algebraic constraints
   *  (\a num_algebraic_cons).
   *
   *  Implementation skeleton:
   *      for (int j=0; j<hdr.num_algebraic_cons; j++) {
   *        AlgConRange bnd;
   *        if (compl_var && compl_var[j]) {
   *          j = compl_var[j]-1;
   *          bnd.k = 0;
   *          if (vlb[j] > negInfinity)
   *            bnd.k = 1;
   *          if (vub[j] < Infinity)
   *            bnd.k |= 2;
   *          assert(bnd.k);
   *          bnd.cvar = j;
   *        } else {
   *          bnd.L = clb[j];
   *          bnd.U = cub[j];
   *        }
   *        cbw.WriteAlgConRange(bnd);
   *      }
   */
  template <class ConBoundsWriter>
  void FeedConBounds(ConBoundsWriter& cbw) {
    for (size_t i=0; i<alg_con_info_.size(); ++i) {          // no constraint permutations
      switch (alg_con_info_[i].GetStaticTypeID()) {
      case StaticItemTypeID::ID_LinConRange: {
        const auto& lcon = *((LinConRange*)(alg_con_info_[i].GetPItem()));
        cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
      } break;
      case StaticItemTypeID::ID_LinConLE: {
        const auto& lcon = *((LinConLE*)(alg_con_info_[i].GetPItem()));
        cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
      } break;
      case StaticItemTypeID::ID_LinConEQ: {
        const auto& lcon = *((LinConEQ*)(alg_con_info_[i].GetPItem()));
        cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
      } break;
      case StaticItemTypeID::ID_LinConGE: {
        const auto& lcon = *((LinConGE*)(alg_con_info_[i].GetPItem()));
        cbw.WriteAlgConRange( AlgConRange{lcon.lb(), lcon.ub()} );
      } break;
      default:
        MP_RAISE("Unknown algebraic constraint type");
      }
    }
  }


  ///////////////////// 6. CONSTRAINTS /////////////////////
  /** Description of constraint \a i
   *    (\a i in 0..num_algebraic_cons+num_logical_cons-1).
   *  With WantNLComments()==true, this is
     *  written to text-format NL as a comment. */
  const char* ConDescription(int ) { return ""; }

  /** Feed the linear part of algebraic constraint \a i.
    * For smooth solvers, should contain entries for all
    * potential nonzeros (Jacobian sparsity pattern).
    *
    *  Implementation skeleton:
    *      if (con_grad[i].size()) {
    *        auto sv = svw.MakeVectorWriter(con_grad[i].size());
    *        for (size_t j=0; j<con_grad.size(); ++j)
    *          sv.Write(con_grad[j].var_index, con_grad[j].coef);
    *      }
    */
  template <class ConLinearExprWriterFactory>
  void FeedLinearConExpr(int i, ConLinearExprWriterFactory& svwf) {
    switch (alg_con_info_[i].GetStaticTypeID()) {
    case StaticItemTypeID::ID_LinConRange: {
      FeedLinearConExpr( *((LinConRange*)(alg_con_info_[i].GetPItem())), svwf);
    } break;
    case StaticItemTypeID::ID_LinConLE: {
      FeedLinearConExpr( *((LinConLE*)(alg_con_info_[i].GetPItem())), svwf);
    } break;
    case StaticItemTypeID::ID_LinConEQ: {
      FeedLinearConExpr( *((LinConEQ*)(alg_con_info_[i].GetPItem())), svwf);
    } break;
    case StaticItemTypeID::ID_LinConGE: {
      FeedLinearConExpr( *((LinConGE*)(alg_con_info_[i].GetPItem())), svwf);
    } break;
    default:
      MP_RAISE("Unknown algebraic constraint type");
    }
  }

  template <class ConLinearExprWriterFactory, class Body, class RhsOrRange>
  void FeedLinearConExpr(
      const AlgebraicConstraint<Body, RhsOrRange>& algcon,
      ConLinearExprWriterFactory& svwf) {
    if (algcon.size()) {
      auto svw = svwf.MakeVectorWriter(algcon.size());
      for (int j=0; j<algcon.size(); ++j)
        svw.Write(GetNewVarIndex( algcon.var(j) ),      // new ordering
                  algcon.coef(j));
    }
  }

  /** Feed nonlinear expression of constraint \a i.
     *  Algebraic constraints (num_algebraic_cons)
     *  come before logical (num_logical_cons).
     *  For linear constraints, the expression should be
   *  constant 0, which is implemented as default.
     */
  template <class ConExprWriter>
  void FeedConExpression(int , ConExprWriter& ew)
  { ew.NPut(0.0); }


  ///////////////////// 7. EXPRESSIONS /////////////////////
  /** Feed native expression.
     *  This method is recursively called from NLWriter,
     *  when Feeder uses ExprWriter::EPut().
     *  Feeder should not call this method itself.
     *
     *  Details of ExprWriter: see NLWriter2.
   */
  template <class ExprWriter>
  void FeedExpr(Expr e, ExprWriter& ) { }


  ///////////////////// 8. PL-SOS CONSTRAINTS ////////////
  /**
     *  The below feature is for AMPL's internal
     *  linearization of piecewise-linear functions.
     *  For user-definable SOS constraints, use suffixes
     *  .sosno/.ref.
     *
     *  The below is a feeder interface
     *  for .sos/.sosref suffixes.
     *  The feeder can provide 3 sparse vectors:
     *  - .sos for variables:
     *    Each nonzero value defines SOS group number.
     *    Negative means SOS Type 2, positive - SOS Type 1.
     *  - .sos for constraints:
     *    Each nonzero value denotes a constraint used in a
     *    linearization of an SOS. The constraint can be deleted
   *    by the solver driver if using solver's SOS.
     *  - .sosref for variables:
     *    SOS weights. Variables participating in an SOS having
     *    zero weights are involved in linearization and can be
     *    deleted if the solver accepts SOS natively.
   *
   *  Implementation:
   *      auto sosv = plsos.StartSOSVars(nvsos);
   *      for (int i=0; i<nvsos; ++i)
   *        sosv.Write(i, vsos[i]);
   *      if (ncsos) {
   *        auto sosc = plsos.StartSOSCons(ncsos);
   *        for ....
   *      }
   *      auto sosrefv = plsos.StartSOSREFVars(ac->nsosref);
   *      ....
    */
  template <class PLSOSWriter>
  void FeedPLSOS(PLSOSWriter& ) { }


  ///////////////////// 9. FUNCTIONS /////////////////////
  /** Function definition. */
  using BaseNLFeeder::FuncDef;

  /** Provide definition
   *  of function \a i, i=0..num_funcs-1. */
  FuncDef Function(int i) { return {}; }


  ///////////////////// 10. RANDOM VARIABLES /////////////////////
  /// Random variables.
  /// Undocumented feature. SNL2006.
  /// Example:
  /// var z >= 0;
  ///	let z.stage := 1;
  ///	var x{0..1, 0..1} random := Uniform(0,2);
  ///	for {i in 0..1, j in 0..1} {let x[i,j].stage := 1;};
  ///	display z.stage, x.stage;
  ///	c: z * sum{i in 0..1, j in 0..1} x[i,j] <= 3 + Sample(Uniform(0,2));
  ///
  /// Feed random variables.
  /// Indexes: num_vars+num_common_exprs
  ///   .. num_vars+num_common_exprs+num_rand_vars-1.
  ///
  /// Implementation skeleton:
  ///     for(j = num_vars+num_common_exprs;
  ///         j < num_vars+num_common_exprs+num_rand_vars; j++) {
  ///       auto ew = rvw.StartRandVar(j, rand_var_comment(j));
  ///       ew.EPut(rand_var_root_expr(j));
  ///     }
  template <class RandVarWriterFactory>
  void FeedRandomVariables(RandVarWriterFactory& ) { }


  ///////////////////// 11. COLUMN SIZES /////////////////////

  /** Jacobian column sizes (including potential nonzeros).
     *  Should feed LP/Jacobian column sizes
     *  for all but the last variable.
     *
     *  This is called before writing Jacobian rows.
   *
   *  Implementation skeleton:
   *      if (WantColumnSizes())
   *        for (int i=0; i < num_vars+num_rand_vars-1; ++i)
   *          csw.Write(col_size[i]);
   */
  template <class ColSizeWriter>
  void FeedColumnSizes(ColSizeWriter& csw) {
    if (WantColumnSizes())
      for (int i=0; i < var_lbs_.size()-1; ++i)        // use old ordering
        csw.Write(mark_data_.col_sizes_orig_[ GetOldVarIndex( i ) ]);
  }


  ///////////////////// 12. INITIAL GUESSES /////////////////////
  /** Initial primal guesses.
   *
   *  Implementation: write all meaningfuls entries (incl. zeros.)
   *      if (ini_guess.size()) {
   *        auto ig = igw.MakeVectorWriter(ini_guess.size());
   *        for (size_t i=0; i<ini_guess.size(); ++i)
   *          ig.Write(ini_guess[i].index_, ini_guess[i].value_);
   *      }
   */
  template <class IGWriter>
  void FeedInitialGuesses(IGWriter& ) { }

  /** Initial dual guesses. */
  template <class IDGWriter>
  void FeedInitialDualGuesses(IDGWriter& ) { }


  ///////////////////// 13. SUFFIXES /////////////////////
  /** Feed suffixes.
     *
     *  For constraints, assume ordering:
     *  first algebraic, then logical.
   *
   *  Implementation: write all non-0 entries (0 is the default.)
   *      while (....) {
   *        auto sw = swf.StartIntSuffix(  // or ...DblSuffix
   *          suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          sw.Write(index[i], value[i]);
   *      }
     */
  template <class SuffixWriterFactory>
  void FeedSuffixes(SuffixWriterFactory& ) { }


  //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
  /** FeedRowAndObjNames:
   *  Provide constraint, then objective names.
   *  Name information is optional.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (i: ....)
   *          wrt << name[i].c_str();
     */
  template <class RowObjNameWriter>
  void FeedRowAndObjNames(RowObjNameWriter& wrt) { }

  /** Provide deleted row names.*/
  template <class DelRowNameWriter>
  void FeedDelRowNames(DelRowNameWriter& ) { }

  /** Provide variable names. */
  template <class ColNameWriter>
  void FeedColNames(ColNameWriter& wrt) {
    if (var_names_.size() && wrt) {
      assert(var_names_.size() == var_lbs_.size());
      for (size_t i=0; i<var_names_.size(); ++i)
        wrt << var_names_[ GetOldVarIndex(i) ];
    }
  }

  /** Provide unused variable names. */
  template <class UnusedVarNameWriter>
  void FeedUnusedVarNames(UnusedVarNameWriter& ) { }

  /** Provide {fixed variable, extra info} pairs.
     *  This includes defined eliminated variables.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (....)
   *          wrt << typename Writer::StrStrValue
     *          { name[i].c_str(), comment[i].c_str() };
     */
  template <class FixedVarNameWriter>
  void FeedFixedVarNames(FixedVarNameWriter& ) { }

  /** Provide {obj name, constant term} pairs.
   *
   *  Implementation:
   *      if (wrt)
   *        for (....)
   *          wrt << typename Writer::StrDblValue
     *          { name[i].c_str(), (double)obj_offset[i] };
     */
  template <class ObjOffsetWriter>
  void FeedObjAdj(ObjOffsetWriter& ) { }


public:
  /// Get new var index for an old var index
  int GetNewVarIndex(int i) const { return mark_data_.var_order_21_[i]; }

  /// Get old var index for a new var index
  int GetOldVarIndex(int i) const { return mark_data_.var_order_12_[i]; }


protected:
  /// For writing NL
  void PrepareModel();

  /// Mark variables.
  /// Now just count binary, integer.
  void MarkVars();

  /// Sort variables
  void SortVars();

  /// Mark and count constraint and variable kinds
  /// (non/linear/integer, range/eqns.)
  void MarkItems();

  NLHeader DoMakeHeader();

  /// Parameters passed when marking variables in an expression tree
  struct ItemMarkingData {
    int n_var_lin_bin_ {0};
    int n_var_lin_int_ {0};

    int n_ranges_ {0};
    int n_eqns_ {0};

    std::size_t n_obj_nz_ {0};
    std::size_t n_con_nz_ {0};

    std::vector< std::pair< int, int > > var_prior_;        // new index -> var weight, orig. index
    std::vector<int> var_order_12_;                         // new index -> old index
    std::vector<int> var_order_21_;                         // old index -> new index

    std::vector<int> col_sizes_orig_;                       // column sizes for original sorting
  };


  /// Parameters passed when submitting a con/obj to the "solver"
  class NLWParams {};


  /// Placeholder for item marker
  template <class Item>
  void MarkItem(const Item& , ItemMarkingData& ) {
    MP_UNSUPPORTED(std::string("MP2NLModelAPI::MarkVars() not supported for ")
                   + typeid(Item).name());
  }

  /// Mark LinConEQ
  void MarkItem(const LinearObjective& lo, ItemMarkingData& prm) {
    mark_data_.n_obj_nz_ += lo.vars().size();        // count nnz
  }

  /// Mark LinConRange
  void MarkItem(const LinConRange& lcr, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += lcr.size();
    Add2ColSizes(lcr.vars());
    if (lcr.lb() > MinusInfinity()
        && lcr.ub() < Infinity()
        && lcr.lb() < lcr.ub())
      ++prm.n_ranges_;
  }

  /// Mark LinConEQ
  void MarkItem(const LinConLE& lcr, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += lcr.size();
    Add2ColSizes(lcr.vars());
  }

  /// Mark LinConEQ
  void MarkItem(const LinConEQ& lcr, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += lcr.size();
    Add2ColSizes(lcr.vars());
    ++prm.n_eqns_;
  }

  /// Mark LinConEQ
  void MarkItem(const LinConGE& lcr, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += lcr.size();
    Add2ColSizes(lcr.vars());
  }

  /// Add to col sizes
  void Add2ColSizes(const std::vector<int>& vars);

  /// Placeholder for bounds writer
  template <class Item>
  void WriteBounds(const Item& , const NLWParams&) {
    MP_UNSUPPORTED(std::string("MP2NLModelAPI::WriteBounds() not supported for ")
                   + typeid(Item).name());
  }


  /// Dispatch con/obj operations: abstract base
  class BasicItemDispatcher {
  public:
    virtual ~BasicItemDispatcher() { }
    /// Construct
    BasicItemDispatcher(MP2NLModelAPI& mapi) : mapi_(mapi) { }
    const MP2NLModelAPI& GetMAPI() const { return mapi_; }
    MP2NLModelAPI& GetMAPI() { return mapi_; }

    virtual void MarkItem(void* pitem, ItemMarkingData& vmp) = 0;
    virtual void WriteBounds(void* pitem, const NLWParams& nlwp) = 0;

  private:
    MP2NLModelAPI& mapi_;
  };


  /// Dispatch con/obj operations: specialization
  template <class Item>
  class ItemDispatcher : public BasicItemDispatcher {
  public:
    /// Construct
    ItemDispatcher(MP2NLModelAPI& mapi)
        : BasicItemDispatcher(mapi) { }
    /// Var marking
    void MarkItem(void* pitem, ItemMarkingData& vmp) override
    { GetMAPI().MarkItem(*(const Item*)pitem, vmp); }
    /// Write bounds
    void WriteBounds(void* pitem, const NLWParams& nlwp) override
    { GetMAPI().WriteBounds(*(const Item*)pitem, nlwp); }
  };


  /// Placeholder for the dispatcher getter
  template <class Item>
  ItemDispatcher<Item>& GetItemDispatcher() { throw 0; }

/// Macro to define an item dispatcher
#define CREATE_ITEM_DISPATCHER(ItemType) \
  ItemDispatcher<ItemType> item_dispatcher_ ## ItemType ## _ { *this }; \
  template <> \
  ItemDispatcher<ItemType>& GetItemDispatcher<ItemType>() \
  { return item_dispatcher_ ## ItemType ## _; }


  /// Item dispatchers
  /// for static item types in the NL format
  CREATE_ITEM_DISPATCHER(LinearObjective)
  CREATE_ITEM_DISPATCHER(NLObjective)

  CREATE_ITEM_DISPATCHER(LinConRange)
  CREATE_ITEM_DISPATCHER(LinConLE)
  CREATE_ITEM_DISPATCHER(LinConEQ)
  CREATE_ITEM_DISPATCHER(LinConGE)
  CREATE_ITEM_DISPATCHER(NLConstraint)
  CREATE_ITEM_DISPATCHER(NLAssignLE)
  CREATE_ITEM_DISPATCHER(NLAssignEQ)
  CREATE_ITEM_DISPATCHER(NLAssignGE)
  CREATE_ITEM_DISPATCHER(NLLogical)
  CREATE_ITEM_DISPATCHER(NLEquivalence)
  CREATE_ITEM_DISPATCHER(NLImpl)
  CREATE_ITEM_DISPATCHER(NLRimpl)
  CREATE_ITEM_DISPATCHER(IndicatorConstraintLinLE)
  CREATE_ITEM_DISPATCHER(IndicatorConstraintLinEQ)
  CREATE_ITEM_DISPATCHER(IndicatorConstraintLinGE)
  CREATE_ITEM_DISPATCHER(SOS1Constraint)
  CREATE_ITEM_DISPATCHER(SOS2Constraint)
  CREATE_ITEM_DISPATCHER(NLComplementarity)


  /// We need item type ID for manual dispatch
  enum class StaticItemTypeID {
    ID_None,
    ID_LinearObjective,
    ID_NLObjective,

    ID_LinConRange,
    ID_LinConLE,
    ID_LinConEQ,
    ID_LinConGE,
    ID_NLConstraint,
    ID_NLAssignLE,
    ID_NLAssignEQ,
    ID_NLAssignGE,
    ID_NLLogical,
    ID_NLEquivalence,
    ID_NLImpl,
    ID_NLRimpl,
    ID_IndicatorConstraintLinLE,
    ID_IndicatorConstraintLinEQ,
    ID_IndicatorConstraintLinGE,
    ID_SOS1Constraint,
    ID_SOS2Constraint,
    ID_NLComplementarity
  };


  /// Expression type IDs for manual dispatch.
  /// These are not NL opcodes:
  /// correspond to our expressions.
  enum class ExpressionTypeID {
    ID_None,

    ID_Lin,
    ID_Quad,

    ID_Abs,
    ID_And,
    ID_Or,
    ID_Exp,
    ID_Log,
    ID_Pow,
    ID_Sin,
    ID_Cos
  };


  /// Constraint/objective/expression info.
  /// We rely on the pointers staying valid.
  class ItemInfo {
  public:
    /// Construct
    ItemInfo (BasicItemDispatcher& disp,
             void* pitem, StaticItemTypeID iid, ExpressionTypeID eid)
        : disp_(disp), p_item_(pitem), itemID_(iid), exprID_(eid) { }
    /// Get dispatcher
    BasicItemDispatcher& GetDispatcher() const { return disp_; }
    /// Get &item
    void* GetPItem() const { return p_item_; }
    /// Get static item type ID, if any
    StaticItemTypeID GetStaticTypeID() const { return itemID_; }
    /// Get expression type ID, if any
    ExpressionTypeID GetExprTypeID() const { return exprID_; }

  private:
    BasicItemDispatcher& disp_;
    void* p_item_;
    StaticItemTypeID itemID_ {StaticItemTypeID::ID_None};
    ExpressionTypeID exprID_ {ExpressionTypeID::ID_None};
  };

  /// Fill static item info
  template <class Item>
  ItemInfo MakeItemInfo(const Item& i, StaticItemTypeID sid) {
    return { GetItemDispatcher<Item>(), (void*)&i, sid, ExpressionTypeID::ID_None };
  }

  /// Fill expression item info
  template <class Item>
  ItemInfo MakeItemInfo(const Item& i, ExpressionTypeID eid) {
    return { GetItemDispatcher<Item>(), (void*)&i, StaticItemTypeID::ID_None, eid };
  }

  /// Create implementations of interfaces of NLSolver
  void CreateInterfaces();

private:
  /// References to the model data.
  /// @note we rely on them staying valid.

  ArrayRef<double> var_lbs_;
  ArrayRef<double> var_ubs_;
  ArrayRef<var::Type> var_types_;
  ArrayRef<const char*> var_names_;                 // can be empty


  /// @todo still need permutations of NL constraints?
  std::vector<ItemInfo> obj_info_, alg_con_info_, log_con_info_, sos_info_;


  ItemMarkingData mark_data_;
  bool hdr_is_current_ {};
  NLHeader hdr_;

  std::unique_ptr<MP2NLSolverIntf> p_nls_;

};

} // namespace mp

#endif // MP2NLMODELAPI_H
