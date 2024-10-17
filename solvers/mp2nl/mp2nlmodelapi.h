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
/// 2. variable expression (represent a (defined) variable),
/// 3. normal expression node.
///
/// @note To create MP2NL_Expr, use MakeMP2NL_... below.
class MP2NL_Expr {
public:
  /// Construct
  MP2NL_Expr(int e=0) : id_(e) { }

	/// Get the expr ID
	int GetID() const { return id_; }

  /// Is this an empty expression?
  bool IsEmptyExpr() const { return !id_; }

  /// Is this expression a variable?
  bool IsVariable() const { return id_ > 0; }

  /// Get variable index
  /// (>=num_vars ==> is a defined variable.)
  int GetVarIndex() const
  { assert(IsVariable()); return id_-1; }

  /// Is a normal expression?
  bool IsExpression() const { return id_ < 0; }

  /// Get expression index
  int GetExprIndex() const
  { assert(IsExpression()); return -id_ - 1; }

private:
	int id_ {0};
};

/// Make an empty expression.
inline MP2NL_Expr MakeEmptyExpr() { return {0}; }

/// Make an expression representing variable \a v.
inline MP2NL_Expr MakeVarExpr(int v) { return {v+1}; }

/// Make ID of a normal expression with index \a i.
inline MP2NL_Expr MakeExprID(int i) { return {-i-1}; }


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
  void InitCustomOptions();

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
	/// Handle expression trees: inherit basic API
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
  /// NLLogical, NLEquivalence, NLImpl, NLRimpl,
  /// and NLObjective.
  ACCEPT_EXPRESSION_INTERFACE(AcceptedButNotRecommended);

  /// Once expressions are supported, need the following
  /// helper methods.
  ///
  /// GetVarExpression(\a i): expression representing variable 0<=i<n_var.
  /// Only called for 'nonlinear' variables.
  Expr GetVarExpression(int i) { return MakeVarExpr(i); }

  /// GetZeroExpr(): constant 0.0 expression.
  /// Can be used to represent empty expression in an NLConstraint.
  Expr GetZeroExpression() { return MakeEmptyExpr(); }

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
  ACCEPT_CONSTRAINT(NLReifEquiv, Recommended, CG_Logical)
  void AddConstraint(const NLReifEquiv& nle);
  /// NL implication: var==1 ==> expression.
  /// This is an expression explicifier in positive context.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLReifImpl, Recommended, CG_Logical)
  void AddConstraint(const NLReifImpl& nle);
  /// NL reverse implication: expression ==> var==1.
  /// This is an expression explicifier in negative context.
  ///
  /// @note Should be 'Recommended'
  ///   whenever logical expressions are accepted.
  /// @note Accessors: GetExpression(nle), GetVariable(nle).
  ACCEPT_CONSTRAINT(NLReifRimpl, Recommended, CG_Logical)
  void AddConstraint(const NLReifRimpl& nle);


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



  /// @brief Accept NLAffineExpr.
  /// Once algebraic expressions are accepted
  /// via NLConstraint, subexpressions might be submitted via
  /// NLAffineExpr and NLQuadExpr.
  ///
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLAffineExpression, Recommended);
  Expr AddExpression(const NLAffineExpression& le);

  /// Accept NLQuadExpr.
  /// @note Use accessors, not methods;
  /// - GetLinSize(le), GetLinCoef(le, i), GetLinTerm(le, i);
  ///   GetQuadSize(le), GetQuadCoef(le, i),
  ///   GetQuadTerm1(le, i), GetQuadTerm2(le, i);
  ///   GetConstTerm(le).
  ACCEPT_EXPRESSION(NLQuadExpression, Recommended);
  Expr AddExpression(const NLQuadExpression& le);

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
  ///
  /// We are not accepting abs() as a constraint, only as expression.
  ACCEPT_CONSTRAINT(AbsConstraint, NotAccepted, CG_General)

  ACCEPT_EXPRESSION(AllDiffExpression, Recommended)
  Expr AddExpression(const AllDiffExpression& absc);

  /// And/Or/Equivalence
  /// @note Use GetNumArguments(expr)
  ///   and GetArgument(expr, i) for i=0..num_args-1.
  ACCEPT_EXPRESSION(AndExpression, Recommended)
  MP2NL_Expr AddExpression(const AndExpression& cc);
  ACCEPT_EXPRESSION(OrExpression, Recommended)
  MP2NL_Expr AddExpression(const OrExpression& dc);
  ACCEPT_EXPRESSION(EquivalenceExpression, Recommended)
  MP2NL_Expr AddExpression(const EquivalenceExpression& dc);

  ACCEPT_EXPRESSION(CondLTExpression, Recommended)
  MP2NL_Expr AddExpression(const CondLTExpression& dc);
  ACCEPT_EXPRESSION(CondLEExpression, Recommended)
  MP2NL_Expr AddExpression(const CondLEExpression& dc);
  ACCEPT_EXPRESSION(CondEQExpression, Recommended)
  MP2NL_Expr AddExpression(const CondEQExpression& dc);
  ACCEPT_EXPRESSION(CondGEExpression, Recommended)
  MP2NL_Expr AddExpression(const CondGEExpression& dc);
  ACCEPT_EXPRESSION(CondGTExpression, Recommended)
  MP2NL_Expr AddExpression(const CondGTExpression& dc);

  ACCEPT_EXPRESSION(IfThenExpression, Recommended)
  MP2NL_Expr AddExpression(const IfThenExpression& dc);
  ACCEPT_EXPRESSION(ImplicationExpression, Recommended)
  MP2NL_Expr AddExpression(const ImplicationExpression& dc);
  ACCEPT_EXPRESSION(NotExpression, Recommended)
  MP2NL_Expr AddExpression(const NotExpression& dc);

  ACCEPT_EXPRESSION(MinExpression, Recommended)
  MP2NL_Expr AddExpression(const MinExpression& dc);
  ACCEPT_EXPRESSION(MaxExpression, Recommended)
  MP2NL_Expr AddExpression(const MaxExpression& dc);


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
  bool WantNLComments() const { return storedOptions_.nl_comments_; }

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
  int ObjType(int i);

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
  void FeedObjGradient(int i, ObjGradWriterFactory& svwf);

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
  void FeedObjExpression(int , ObjExprWriter& ew);


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
  void FeedVarBounds(VarBoundsWriter& vbw);


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
  void FeedConBounds(ConBoundsWriter& cbw);

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
  void FeedLinearConExpr(int i, ConLinearExprWriterFactory& svwf);

  template <class ExprBody, class ConLinearExprWriterFactory>
  void FeedLinearConBody(
      const ExprBody& algcon,
      ConLinearExprWriterFactory& svwf);

  /** Feed nonlinear expression of constraint \a i.
     *  Algebraic constraints (num_algebraic_cons)
     *  come before logical (num_logical_cons).
     *  For linear constraints, the expression should be
   *  constant 0, which is implemented as default.
     */
  template <class ConExprWriter>
  void FeedConExpression(int , ConExprWriter& ew);

  /// Feed alg con expr
  template <class ConExprWriter>
  void FeedAlgConExpression(int , ConExprWriter& ew);

  /// Feed logical con expr
  template <class ConExprWriter>
  void FeedLogicalConExpression(int , ConExprWriter& ew);


  ///////////////////// 7. EXPRESSIONS /////////////////////
  /** Feed native expression.
     *  This method is recursively called from NLWriter,
     *  when Feeder uses ExprWriter::EPut().
     *  Feeder should not call this method
     *  to write subtrees below the root expression.
     *
     *  Details of ExprWriter: see NLWriter2.
   */
  template <class ExprWriter>
  void FeedExpr(Expr e, ExprWriter& );

  /// Write opcode(s)
  template <class ExprWriter>
  void FeedOpcode(Expr e, ExprWriter& );

  /// Feed NLAffineExpression or NLQuadExpression
  template <class AlgMPExpr, class ExprWriter>
  void FeedAlgebraic(const AlgMPExpr& e, ExprWriter& ew);

  /// Feed comparison
  template <class CondMPExpr, class ExprWriter>
  void FeedRelational(const CondMPExpr& e, ExprWriter ew);

  /// Feed NLLogical
  template <class ExprWriter>
  void FeedNLLogical(const NLLogical& e, ExprWriter& ew);

  /// Feed NLBaseReif
  template <int sense, class ExprWriter>
  void FeedReification(const NLBaseReif<sense>& e, ExprWriter& ew);

  /// Write opcode arguments.
  /// Parameters written after the expression arguments.
  /// @param aw: argument writer produced by OPutN() or similar.
  template <class MPExpr, class ArgWriter>
  void FdArgs(const MPExpr& e, ArgWriter aw);

  /// Write opcode logical arguments.
  /// Variables are written as "var==1".
  /// Parameters written after the expression arguments.
  /// @param aw: argument writer produced by OPutN() or similar.
  template <class MPExpr, class ArgWriter>
  void FdLogicArgs(const MPExpr& e, ArgWriter aw);

  /// Write IfThen's arguments.
  /// @todo the then/else arguments might need to be logical...
  template <class MPExpr, class ArgWriter>
  void FdIfArgs(const MPExpr& e, ArgWriter aw);

  /// Write a logical argument.
  /// Variables are written as "var==1".
  template <class ArgWriter>
  void FeedLogicalExpression(MP2NL_Expr mp2nle, ArgWriter& aw);


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
  void FeedColumnSizes(ColSizeWriter& csw);


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
  void FeedInitialGuesses(IGWriter& );

  /** Initial dual guesses. */
  template <class IDGWriter>
  void FeedInitialDualGuesses(IDGWriter& );


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
  /// @todo SOS constraints
  template <class SuffixWriterFactory>
  void FeedSuffixes(SuffixWriterFactory& );


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
  void FeedRowAndObjNames(RowObjNameWriter& wrt);

  /** Provide deleted row names.*/
  template <class DelRowNameWriter>
  void FeedDelRowNames(DelRowNameWriter& ) { }

  /** Provide variable names. */
  template <class ColNameWriter>
  void FeedColNames(ColNameWriter& wrt);

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

  /// Item name
  template <class Item>
  const char* GetItemName(const Item& item) const
  { return item.name(); }

  /// Item name
  template <class FuncCon>
  const char* GetItemName(const ExprWrapper<FuncCon>& item) const
  { return item.GetFlatConstraint().name(); }

  /// NL/SOL file stub
  const std::string& GetFileStub() const
  { return storedOptions_.stub_; }

protected:
  /// For writing NL
  void PrepareModel();

  /// Map expressions.
  void MapExpressions();

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
    MP_UNSUPPORTED(std::string("MP2NLModelAPI::MarkItem() not supported for ")
                   + typeid(Item).name());
  }

  /// Mark linear part of any objective
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

  /// Mark NLconstraint
  void MarkItem(const NLConstraint& nlc, ItemMarkingData& prm) {
    MarkItem(nlc.GetMainCon(), prm);
  }

  /// Mark NLAssign
  void MarkItem(const NLAssignLE& nlc, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += 1;
    Add2ColSizes({ nlc.GetVar() });
  }
  /// Mark NLAssign
  void MarkItem(const NLAssignEQ& nlc, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += 1;
    Add2ColSizes({ nlc.GetVar() });
    ++prm.n_eqns_;
  }
  /// Mark NLAssign
  void MarkItem(const NLAssignGE& nlc, ItemMarkingData& prm) {
    mark_data_.n_con_nz_ += 1;
    Add2ColSizes({ nlc.GetVar() });
  }

  /// Add to col sizes
  void Add2ColSizes(const std::vector<int>& vars);

  /// Placeholder for bounds writer
  template <class Item>
  void WriteBounds(const Item& , const NLWParams&) {
    MP_UNSUPPORTED(std::string("MP2NLModelAPI::WriteBounds() not supported for ")
                   + typeid(Item).name());
  }


  /// We need item type ID for manual dispatch.
  /// Manual dispatch is necessary to pass new types,
  /// e.g., from NLWriter.
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
    ID_NLReifEquiv,
    ID_NLReifImpl,
    ID_NLReifRimpl,
    ID_IndicatorConstraintLinLE,
    ID_IndicatorConstraintLinEQ,
    ID_IndicatorConstraintLinGE,
    ID_SOS1Constraint,
    ID_SOS2Constraint,
    ID_NLComplementarity
  };


  /// Expression type IDs for manual dispatch.
  /// These are not NL opcodes:
  /// instead, they correspond to our expressions.
  enum class ExpressionTypeID {
    ID_None,

    ID_NLAffine,
    ID_NLQuad,

    ID_AllDiff,
    ID_And,
    ID_Not,
    ID_Or,
    ID_Equivalence,
    ID_IfThen,
    ID_Implication,

    ID_CondLT,
    ID_CondLE,
    ID_CondEQ,
    ID_CondGE,
    ID_CondGT,

    ID_Abs,
    ID_Min,
    ID_Max,

    ID_Exp,
    ID_Log,
    ID_Pow,
    ID_Sin,
    ID_Cos
  };


  /// Normally dispatch con/obj operations: abstract base.
  /// This kind of dispatch only works
  /// with known types.
  /// To pass new types, we'll use manual dispatch.
  class BasicItemDispatcher {
  public:
    virtual ~BasicItemDispatcher() { }

    /// Construct
    BasicItemDispatcher(MP2NLModelAPI& mapi) : mapi_(mapi) { }
    const MP2NLModelAPI& GetMAPI() const { return mapi_; }
    MP2NLModelAPI& GetMAPI() { return mapi_; }

    /// ItemName
    virtual const char* GetName(void* pitem) = 0;

    virtual void MarkItem(void* pitem, ItemMarkingData& vmp) = 0;
    virtual void WriteBounds(void* pitem, const NLWParams& nlwp) = 0;
    virtual void MarkExprTree(void* pitem) = 0;

    /// Placeholder for the item category getter:
    /// static item vs expression
    virtual bool IsItemTypeStatic() = 0;

    /// Placeholder for the StaticItemTypeID getter
    virtual StaticItemTypeID GetStaticItemTypeID() = 0;

    /// Placeholder for the ExpressionTypeID getter
    virtual ExpressionTypeID GetExpressionTypeID() = 0;

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

    /// ItemName
    const char* GetName(void* pitem) override
    { return GetMAPI().GetItemName(*(const Item*)pitem); }

    /// Var marking
    void MarkItem(void* pitem, ItemMarkingData& vmp) override
    { GetMAPI().MarkItem(*(const Item*)pitem, vmp); }

    /// Write bounds
    void WriteBounds(void* pitem, const NLWParams& nlwp) override
    { GetMAPI().WriteBounds(*(const Item*)pitem, nlwp); }

    /// Add + mark expressions
    void MarkExprTree(void* pitem) override
    { GetMAPI().MapExprTreeFromItem(*(const Item*)pitem); }

    /// Item category getter:
    /// static item vs expression
    bool IsItemTypeStatic() override
    { return MP2NLModelAPI::IsItemTypeStatic<Item>(); }

    /// Placeholder for the StaticItemTypeID getter
    StaticItemTypeID GetStaticItemTypeID() override
    { return MP2NLModelAPI::GetStaticItemTypeID<Item>(); }

    /// Placeholder for the ExpressionTypeID getter
    ExpressionTypeID GetExpressionTypeID() override
    { return MP2NLModelAPI::GetExpressionTypeID<Item>(); }
  };


  /// Placeholder for the dispatcher getter
  template <class Item>
  ItemDispatcher<Item>& GetItemDispatcher() { throw 0; }

  /// Placeholder for the item category getter:
  /// static item vs expression
  template <class Item>
  static bool IsItemTypeStatic() { throw 0; }

  /// Placeholder for the StaticItemTypeID getter
  template <class Item>
  static StaticItemTypeID GetStaticItemTypeID() { throw 0; }

  /// Placeholder for the ExpressionTypeID getter
  template <class Item>
  static ExpressionTypeID GetExpressionTypeID() { throw 0; }


/// Macro to define an item dispatcher
#define CREATE_STATIC_ITEM_DISPATCHER(ItemType) \
  ItemDispatcher<ItemType> item_dispatcher_ ## ItemType ## _ { *this }; \
  template <> \
  ItemDispatcher<ItemType>& GetItemDispatcher<ItemType>() \
  { return item_dispatcher_ ## ItemType ## _; } \
  template <> \
  bool IsItemTypeStatic<ItemType>() { return true; } \
  template <> \
  StaticItemTypeID GetStaticItemTypeID<ItemType>() \
  { return StaticItemTypeID::ID_ ## ItemType; }


  /// Item dispatchers
  /// for static item types in the NL format
  CREATE_STATIC_ITEM_DISPATCHER(LinearObjective)
  CREATE_STATIC_ITEM_DISPATCHER(NLObjective)

  CREATE_STATIC_ITEM_DISPATCHER(LinConRange)
  CREATE_STATIC_ITEM_DISPATCHER(LinConLE)
  CREATE_STATIC_ITEM_DISPATCHER(LinConEQ)
  CREATE_STATIC_ITEM_DISPATCHER(LinConGE)
  CREATE_STATIC_ITEM_DISPATCHER(NLConstraint)
  CREATE_STATIC_ITEM_DISPATCHER(NLAssignLE)
  CREATE_STATIC_ITEM_DISPATCHER(NLAssignEQ)
  CREATE_STATIC_ITEM_DISPATCHER(NLAssignGE)
  CREATE_STATIC_ITEM_DISPATCHER(NLLogical)
  CREATE_STATIC_ITEM_DISPATCHER(NLReifEquiv)
  CREATE_STATIC_ITEM_DISPATCHER(NLReifImpl)
  CREATE_STATIC_ITEM_DISPATCHER(NLReifRimpl)
  CREATE_STATIC_ITEM_DISPATCHER(IndicatorConstraintLinLE)
  CREATE_STATIC_ITEM_DISPATCHER(IndicatorConstraintLinEQ)
  CREATE_STATIC_ITEM_DISPATCHER(IndicatorConstraintLinGE)
  CREATE_STATIC_ITEM_DISPATCHER(SOS1Constraint)
  CREATE_STATIC_ITEM_DISPATCHER(SOS2Constraint)
  CREATE_STATIC_ITEM_DISPATCHER(NLComplementarity)


/// Macro to define an expression dispatcher
#define CREATE_EXPRESSION_DISPATCHER(ItemType) \
  ItemDispatcher<ItemType ## Expression> \
    item_dispatcher_ ## ItemType ## _ { *this }; \
  template <> \
  ItemDispatcher<ItemType ## Expression>& \
  GetItemDispatcher<ItemType ## Expression>() \
  { return item_dispatcher_ ## ItemType ## _; } \
  template <> \
  bool IsItemTypeStatic<ItemType ## Expression>() \
  { return false; } \
  template <> \
  ExpressionTypeID GetExpressionTypeID<ItemType ## Expression>() \
  { return ExpressionTypeID::ID_ ## ItemType; }


  CREATE_EXPRESSION_DISPATCHER(NLAffine)
  CREATE_EXPRESSION_DISPATCHER(NLQuad)

  CREATE_EXPRESSION_DISPATCHER(And)
  CREATE_EXPRESSION_DISPATCHER(Or)
  CREATE_EXPRESSION_DISPATCHER(Equivalence)

  CREATE_EXPRESSION_DISPATCHER(AllDiff)
  CREATE_EXPRESSION_DISPATCHER(CondLT)
  CREATE_EXPRESSION_DISPATCHER(CondLE)
  CREATE_EXPRESSION_DISPATCHER(CondEQ)
  CREATE_EXPRESSION_DISPATCHER(CondGE)
  CREATE_EXPRESSION_DISPATCHER(CondGT)
  CREATE_EXPRESSION_DISPATCHER(IfThen)
  CREATE_EXPRESSION_DISPATCHER(Implication)
  CREATE_EXPRESSION_DISPATCHER(Not)

  CREATE_EXPRESSION_DISPATCHER(Abs)
  CREATE_EXPRESSION_DISPATCHER(Min)
  CREATE_EXPRESSION_DISPATCHER(Max)

  CREATE_EXPRESSION_DISPATCHER(Exp)
  CREATE_EXPRESSION_DISPATCHER(Log)
  CREATE_EXPRESSION_DISPATCHER(Pow)
  CREATE_EXPRESSION_DISPATCHER(Sin)
  CREATE_EXPRESSION_DISPATCHER(Cos)


  /// Constraint/objective/expression info.
  /// We rely on the pointers staying valid.
  class ItemInfo {
  public:
    /// Construct
    ItemInfo (BasicItemDispatcher& disp, void* pitem,
             bool fLogical
             //, StaticItemTypeID iid, ExpressionTypeID eid
             )
        : disp_(disp), p_item_(pitem), f_logical_(fLogical)
        // no storing, itemID_(iid), exprID_(eid)
    { }
    /// Get dispatcher
    BasicItemDispatcher& GetDispatcher() const { return disp_; }
    /// Get &item
    void* GetPItem() const { return p_item_; }
    /// Is logical?
    bool IsLogical() const { return f_logical_; }
    /// Is a static type?
    bool IsItemTypeStatic() const
    { return GetDispatcher().IsItemTypeStatic(); }
    /// Get static item type ID, if any
    StaticItemTypeID GetStaticTypeID() const
    { return GetDispatcher().GetStaticItemTypeID(); }
    /// Get expression type ID, if any
    ExpressionTypeID GetExprTypeID() const
    { return GetDispatcher().GetExpressionTypeID(); }

  private:
    BasicItemDispatcher& disp_;
    void* p_item_;
    bool f_logical_ {};
    // StaticItemTypeID itemID_ {StaticItemTypeID::ID_None};
    // ExpressionTypeID exprID_ {ExpressionTypeID::ID_None};
  };

  /// Fill static item info
  template <class Item>
  ItemInfo MakeItemInfo(
      const Item& i, StaticItemTypeID , bool fLogical) {
    return { GetItemDispatcher<Item>(), (void*)&i, fLogical
            // , sid, ExpressionTypeID::ID_None
    };
  }

  /// Fill expression item info
  template <class Item>
  ItemInfo MakeItemInfo(
      const Item& i, ExpressionTypeID , bool fLogical) {
    return { GetItemDispatcher<Item>(), (void*)&i, fLogical
            //, StaticItemTypeID::ID_None, eid
    };
  }

  /// Manually dispatch a static item
  /// to a lambda(const auto& item).
  template <class Lambda>
  void DispatchStaticItem(const ItemInfo& info, Lambda lambda);

  /// Map expressions from a single item info
  /// @note needs to be in .h to be instantiated/inlined,
  ///   at least for Clang 16.
  void MapExprTreeFromItemInfo(const ItemInfo& info)
  { info.GetDispatcher().MarkExprTree(info.GetPItem()); }

  /// Map expressions from a single item (type-safe)
  /// @note needs to be in .h to be instantiated/inlined,
  ///   at least for Clang 16.
  template <class Item>
  void MapExprTreeFromItem(const Item& item) {
    auto mp2nlexpr
        = GetExpression(item);   // so that AddExpression() is called
    RegisterExpression(mp2nlexpr);     // mark the tree root
  }

  /// Reuse GetExpression()
  using BaseModelAPI::GetExpression;

  /// Supply empty expression for other types
  /// @note needs to be in .h to be instantiated/inlined,
  ///   at least for Clang 16.
  template <class Item>
  MP2NL_Expr GetExpression(const Item& )
  { return GetZeroExpression(); }

  /// Create implementations of interfaces of NLSolver
  void CreateInterfaces();


protected:
  /// Every new expression when adding.
  /// Just counts them.
  void RegisterExpression(MP2NL_Expr expr);

  /// Count expression depending on its kind.
  void CountExpression(MP2NL_Expr expr);

  /// Single template code to add any expression
  /// where eid is provided
  template <class Expr>
  MP2NL_Expr AddExpression(
      const Expr& expr, ExpressionTypeID eid);


  /// Make, store and return MP2NL_Expr for a given expession
  template <class Expr>
  MP2NL_Expr StoreMP2NLExprID(
      const Expr& expr, ExpressionTypeID eid);

  /// Variable name
  const char* GetVarName(int v) const {
    assert(v>=0 && v<(int)var_lbs_.size());
    return
        (int)var_names_.size() > v ? var_names_[v] : nullptr;
  }


private:
  /// References to the model data.
  /// @note we rely on them staying valid.

  ArrayRef<double> var_lbs_;
  ArrayRef<double> var_ubs_;
  ArrayRef<var::Type> var_types_;
  ArrayRef<const char*> var_names_;                 // can be empty


  /// @todo still need permutations of NL constraints?
  std::vector<ItemInfo> obj_info_, alg_con_info_, log_con_info_, sos_info_;

  std::vector<ItemInfo> expr_info_;
  std::vector<int>   expr_counter_;   // usage counter

  ItemMarkingData mark_data_;
  bool hdr_is_current_ {};
  NLHeader hdr_;

  std::unique_ptr<MP2NLSolverIntf> p_nls_;

  /// These options are stored in the class
  struct Options {
    std::string stub_;
    int nl_comments_ {};
    int nl_format_text_ {};
  };
  Options storedOptions_;
};

} // namespace mp

#endif // MP2NLMODELAPI_H
