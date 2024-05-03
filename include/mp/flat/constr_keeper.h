#ifndef CONSTRAINT_KEEPER_H
#define CONSTRAINT_KEEPER_H

#include <deque>
#include <unordered_map>
#include <functional>
#include <cmath>

#include "mp/common.h"
#include "mp/format.h"
#include "mp/env.h"
#include "mp/utils-math.h"
#include "mp/utils-file.h"
#include "mp/util-json-write.hpp"

#include "mp/flat/model_api_base.h"
#include "mp/flat/constr_hash.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_eval.h"

namespace mp {


/// Converters handling custom constraints should derive from
class BasicFlatConverter;


/// Interface for an array of constraints of certain type
class BasicConstraintKeeper {
public:
  /// Destructor
  virtual ~BasicConstraintKeeper() { }

  /// Constructor
  BasicConstraintKeeper(
      pre::BasicValuePresolver& pres,
      const char* nm, const char* optN) :
    value_node_(pres, nm),
    constr_name_(nm), solver_opt_nm_(optN) { }

  /// Constraint type
  using ConstraintType = BasicConstraint;

  /// Constraint keeper description
  virtual const std::string& GetDescription() const = 0;


  /// Get context of contraint \a i
  virtual Context GetContext(int i) const = 0;

  /// Propagate expression result of constraint \a i top-down
  virtual void PropagateResult(BasicFlatConverter& cvt,
                               int i,
                               double lb, double ub, Context ctx) = 0;

  /// Result variable of constraint \a i. Returns -1 if none
  virtual int GetResultVar(int i) const = 0;

  /// Convert all new items of this constraint.
  /// This normally dispatches conversion (decomposition)
  ///  to the Converter
  /// @return whether any converted
  virtual bool ConvertAllNewWith(BasicFlatConverter& cvt) = 0;

  /// Mark whether we could keep result vars
  virtual void MarkExprResultVars(BasicFlatConverter& cvt) = 0;

  /// Then, mark flat constraint arguments as vars
  virtual void MarkArguments(BasicFlatConverter& cvt) = 0;

  /// Convert to use expressions
  virtual void ConvertWithExpressions(BasicFlatConverter& cvt) = 0;

  /// Query (user-chosen) acceptance level.
  /// This is "combined" for constraint or expression
  ConstraintAcceptanceLevel GetChosenAcceptanceLevel() const {
    if (acceptance_level_<0) {      // not initialized
      std::array<int, 5> alv = {0, 1, 2, 1, 2};
      acceptance_level_ = alv.at(acc_level_item_);
    }
    return ConstraintAcceptanceLevel(acceptance_level_);
  }

  /// Query (user-chosen) expression acceptance level.
  ExpressionAcceptanceLevel GetChosenAcceptanceLevelEXPR() const {
    if (acc_level_expr_<0) {      // not initialized
      std::array<int, 5> alv = {0, 0, 0, 1, 2};
      acc_level_expr_ = alv.at(acc_level_item_);
    }
    return ExpressionAcceptanceLevel(acc_level_expr_);
  }

  /// Converter's ability to convert the constraint type
  virtual bool IfConverterConverts(
      BasicFlatConverter& cvt ) const = 0;

  /// ModelAPI's acceptance level for the constraint type.
  /// This should not be used directly, instead:
  /// GetChosenAcceptanceLevel()
  virtual ConstraintAcceptanceLevel GetModelAPIAcceptance(
      const BasicFlatModelAPI& ) const = 0;

  /// ModelAPI's acceptance level for the expression type.
  /// This should not be used directly, instead:
  /// GetChosenAcceptanceLevelEXPR()
  virtual ExpressionAcceptanceLevel GetModelAPIAcceptanceEXPR(
      const BasicFlatModelAPI& ) const = 0;

  /// Acceptance level of the overall expression interface in the ModelAPI
  virtual ExpressionAcceptanceLevel GetModelAPIAcceptance_EXPR_INTF(
      const BasicFlatModelAPI& ba) const = 0;

  /// Constraint type_info
  virtual const std::type_info& GetTypeInfo() const =0;

  /// Backend's group number for the constraint type
  virtual int GetConstraintGroup(const BasicFlatModelAPI& ) const = 0;

  /// Report how many will be added to Backend
  virtual int GetNumberOfAddable() const = 0;

  /// This adds all unbridged items to the backend (without conversion)
  virtual void AddUnbridgedToBackend(
      BasicFlatModelAPI& be, const std::vector<std::string>* vnames) = 0;

  /// This logs the constraint group
  virtual void LogConstraintGroup(BasicFlatModelAPI& be) = 0;

  /// Value presolve node, const
  const pre::ValueNode& GetValueNode() const { return value_node_; }

  /// Value presolve node
  pre::ValueNode& GetValueNode() { return value_node_; }

  /// Create ValueNode range pointer: add n elements
  pre::NodeRange AddValueNodeRange(int n=1)
  { return GetValueNode().Add(n); }

  /// Create ValueNode range pointer: select n elements at certain pos
  pre::NodeRange SelectValueNodeRange(int pos, int n=1)
  { return GetValueNode().Select(pos, n); }

  /// Constraint type name, e.g., 'AbsConstraint'
  const char* GetConstraintName() const { return constr_name_; }

  /// Expression type, or, if appropriate, constraint type name,
  /// e.g., 'Abs'
  virtual const char* GetExprOrConstraintName() const =0;

  /// Acceptance option names
  virtual const char* GetAcceptanceOptionNames() const
  { return solver_opt_nm_; }

  /// Constraint type short name.
  /// Ideally should be in the constraint itself,
  /// but currently we derive it from acceptance options.
  virtual const char* GetShortTypeName() const;

  /// See what options are available for this constraint:
  /// whether it is accepted natively by ModelAPI,
  /// as flat constraint or expression.
  /// Add acceptance option(s) "acc:...".
  /// Populate constraint list for -c output.
  /// @note This should be called before using the class.
  void ConsiderAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env) {
    DoAddAcceptanceOptions(cvt, ma, env);
    DoPopulateConstraintList(cvt, ma, env);  // for -c option
  }

  /// Mark as bridged. Use index only.
  virtual void MarkAsBridged(int i) = 0;

  /// Mark as unused. Use index only.
  virtual void MarkAsUnused(int i) = 0;

  /// Is constraint \a i unused?
  virtual bool IsUnused(int i) const = 0;

  /// Copy names from ValueNodes
  virtual void CopyNamesFromValueNodes() = 0;

  /// Compute result for constraint \a i
  /// (for functional constraints)
  virtual double ComputeValue(int i, const VarInfoRecomp& ) = 0;

  /// Compute violations
  virtual void ComputeViolations(SolCheck& ) = 0;

  /// Set logger
  void SetLogger(BasicLogger* lg) { exporter_=lg; }
  /// Get logger, if provided and open and ok.
  BasicLogger* GetLogger() const {
    return exporter_ && exporter_->IsOpen()
        ? exporter_ : nullptr;
  }

protected:
  void DoAddAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env);
  /// For -c
  void DoPopulateConstraintList(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env);


private:
  pre::ValueNode value_node_;
  const char* const constr_name_;
  const char* const solver_opt_nm_;
  mutable std::string type_name_short_;
  mutable int acceptance_level_ {-1};     // combined, for either con or expr
  int acc_level_item_ {0};                // item, corresp. to the solver option 0..4
  mutable int acc_level_expr_ {-1};       // expression only
  BasicLogger* exporter_{};
};


/// Full id of a constraint: CK + index
/// This helper class is parameterized by the Keeper
template <class ConstraintKeeper>
struct ConstraintLocationHelper {
  /// Default constructor: no valid Id
  ConstraintLocationHelper() = default;
  /// Normal constructor
  ConstraintLocationHelper(ConstraintKeeper* pck, int i) noexcept :
    pck_(pck), index_(i) { }

  /// Checks if we store a constraint's location
  operator bool() const { return HasId(); }
  /// Checks if we store a constraint's location
  bool HasId() const { return nullptr!=pck_; }
  /// High-level getter
  int GetResultVar() const { return GetCK()->GetResultVar(GetIndex()); }
  /// High-level getter
  const typename ConstraintKeeper::ConstraintType&
  GetConstraint() const { return GetCK()->GetConstraint(GetIndex()); }
  /// High-level getter
  typename ConstraintKeeper::ConstraintType&
  GetConstraint() { return GetCK()->GetConstraint(GetIndex()); }

  /// Get Keeper
  ConstraintKeeper* GetCK() const { assert(HasId()); return pck_; }
  /// Get index
  int GetIndex() const { return index_; }

  /// Set Keeper
  void SetCK(ConstraintKeeper* pck) { pck_ = pck; }
  /// Set index
  void SetIndex(int i) { index_ = i; }

  ConstraintKeeper* pck_ = nullptr;
  int index_ = 0;        // constraint index
};


/// Without constraint type
using AbstractConstraintLocation =
  ConstraintLocationHelper<BasicConstraintKeeper>;


/// Converters handling custom constraints should derive from
class BasicFlatConverter {
public:
  /// Default conversion priority
  static constexpr double ConstraintCvtPriority(BasicConstraint*) { return 1.0; }

  /// Derived converter classes have to tell C++ to use
  /// default handlers if they need them
  /// when they overload Convert() etc, due to C++ name hiding
#define USE_BASE_CONSTRAINT_CONVERTERS(BaseConverter) \
  using BaseConverter::PreprocessConstraint; \
  using BaseConverter::PropagateResult; \
  using BaseConverter::IfHasCvt_impl; \
  using BaseConverter::IfNeedsCvt_impl; \
  using BaseConverter::Convert


  /// For Common Subexpression Elimination, we can use maps
  /// This stub returns empty Id
  int MapFind(const BasicConstraint& ) { return -1; }

  /// Returns false when we do have a map and entry duplicated
  /// (should not happen).
  /// Can be conveniently overloaded
  template <class Con>
  bool MapInsert(const Con& , int ) { return true; }

  /// Similarly to Convert(),
  /// need to 'using' base class' map accessors in the Converter
#define USE_BASE_MAP_FINDERS(BaseConverter) \
  using BaseConverter::MapFind; \
  using BaseConverter::MapInsert; \
  template <class Constraint> \
  using ConstraintLocation = \
    ConstraintLocationHelper< \
      ConstraintKeeper< Impl, ModelAPI, Constraint > >;


  /// Value of Pi
  static constexpr double Pi() { return 3.14159265358979; }

  /// Infinity
  static constexpr double Infty() { return INFINITY; }
  /// -Infinity
  static constexpr double MinusInfty() { return -INFINITY; }
  /// Pract inf
	static constexpr double PracticallyInf() { return 1e20; }
  /// Pract -inf
	static constexpr double PracticallyMinusInf() { return -1e20; }
};


/// Specialize ConstraintKeeper for a given constraint type
/// to store an array of such constraints
template <class Converter, class Backend, class Constraint>
class ConstraintKeeper final
    : public BasicConstraintKeeper {
public:
  /// Constructor, adds this CK to the provided ConstraintManager
  /// Requires the CM to be already constructed
  ConstraintKeeper(Converter& cvt, const char* nm, const char* optnm) :
      BasicConstraintKeeper(cvt.GetValuePresolver(), nm, optnm), cvt_(cvt)
  {
    GetValueNode().SetName(GetShortTypeName());  // change value node name
    GetConverter().AddConstraintKeeper(*this, ConversionPriority());
  }

  /// Constraint type
  using ConstraintType = Constraint;

  /// The corresponding flat expression type
  using FlatExprType = ExprWrapper<Constraint>;

  /// Expression type, or, if appropriate, constraint type name,
  /// e.g., 'Abs'
  const char* GetExprOrConstraintName() const override
  { return Constraint::GetTypeName(); }

  /// Constrint Keeper description
  const std::string& GetDescription() const override
  { return desc_; }

  /// Assume Converter has the Backend
  Backend& GetBackend(BasicFlatConverter& cvt)
  { return static_cast<Converter&>(cvt).GetModelAPI(); }

  /// Add a pre-constructed constraint (or just arguments)
  /// @return index of the new constraint
  template <class... Args>
  int AddConstraint(int d, Args&&... args)
  {
    cons_.emplace_back( d, std::move(args)... );
    ExportConstraint(cons_.size()-1, cons_.back());
    return cons_.size()-1;
  }

  /// Get const constraint \a i
  const Constraint& GetConstraint(int i) const
  { assert(check_index(i)); return cons_[i].GetCon(); }

  /// Get constraint \a i
  Constraint& GetConstraint(int i)
  { assert(check_index(i)); return cons_[i].GetCon(); }

  /// Get constraint depth in the reformulation tree
  int GetConstraintDepth(int i) const
  { assert(check_index(i)); return cons_[i].GetDepth(); }

  /// Get context of contraint \a i
  Context GetContext(int i) const override
  { assert(check_index(i)); return cons_[i].GetCon().GetContext(); }

  /// Propagate expression result of constraint \a i top-down
  void PropagateResult(BasicFlatConverter& cvt,
                       int i,
                       double lb, double ub, Context ctx) override {
    try {
      static_cast<Converter&>(cvt).PropagateResult(
            GetConstraint(i), lb, ub, ctx);
    } catch (const std::exception& exc) {
      MP_RAISE(Converter::GetTypeName() +
               std::string(": propagating result for constraint ") +
               std::to_string(i) + " of type '" +
               Constraint::GetTypeName() +
               "':  " + exc.what());
    }
  }

  /// Result variable of constraint \a i. Returns -1 if none
  int GetResultVar(int i) const override
  { assert(check_index(i)); return cons_[i].GetCon().GetResultVar(); }

  /// Conversion priority. Uses that from Converter
  double ConversionPriority() const
  { return Converter::ConstraintCvtPriority((Constraint*)nullptr); }

  /// Convert all new items of this constraint.
  /// This normally dispatches conversion (decomposition) to the Converter
  /// @return whether any converted
  bool ConvertAllNewWith(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    MP_UNUSED(cvt);
    try {
      return ConvertAllFrom(i_cvt_last_);
    } catch (const std::exception& exc) {
      MP_RAISE(Converter::GetTypeName() + std::string(": ")
                             + exc.what());
    }
    return false;
  }

  /// Mark whether we could result vars of functional constraints
  /// as vars, vs using these constraints as expressions
  void MarkExprResultVars(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    DoMarkForResultVars();
  }

  /// Then, mark arguments of flat constraints as proper vars
  void MarkArguments(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    DoMarkForArguments();
  }

  /// Convert to use expressions
  void ConvertWithExpressions(BasicFlatConverter& cvt) override {
    assert(&cvt == &GetConverter());         // Using the same Converter
    DoCvtWithExprs();
  }

  /// Converter's ability to convert the constraint type
  bool IfConverterConverts(
      BasicFlatConverter& cvt ) const override {
    return static_cast<Converter&>(cvt).
        IfHasConversion((const Constraint*)nullptr);
  }

  /// Acceptance level of this constraint type in the ModelAPI
  ConstraintAcceptanceLevel GetModelAPIAcceptance(
      const BasicFlatModelAPI& ba) const override {
    return
        static_cast<const Backend&>( ba ).
        AcceptanceLevel((Constraint*)nullptr);
  }

  /// Acceptance level of the corresponding expression type in the ModelAPI
  ExpressionAcceptanceLevel GetModelAPIAcceptanceEXPR(
      const BasicFlatModelAPI& ba) const override {
    return
        static_cast<const Backend&>( ba ).
        AcceptanceLevel((FlatExprType*)nullptr);
  }

  /// Acceptance level of the overall expression interface in the ModelAPI
  ExpressionAcceptanceLevel GetModelAPIAcceptance_EXPR_INTF(
      const BasicFlatModelAPI& ba) const override {
    return
        static_cast<const Backend&>( ba ).
        ExpressionInterfaceAcceptanceLevel();
  }

  /// Constraint type_info
  const std::type_info& GetTypeInfo() const override
  { return typeid(ConstraintType); }

  /// Report how many will be added to Backend
  int GetNumberOfAddable() const override {
    return (int)cons_.size()-n_bridged_or_unused_;
  }

  /// Group number of this constraint type in the Backend.
  /// This is needed for pre- / postsolve to group solution values
  int GetConstraintGroup(const BasicFlatModelAPI& ba) const override {
    return static_cast<const Backend&>( ba ).
        GroupNumber((Constraint*)nullptr);
  }

  /// Add remaining constraints to Backend
  void AddUnbridgedToBackend(
      BasicFlatModelAPI& be,
      const std::vector<std::string>* pvnam) override {
    try {
      AddAllUnbridged(be, pvnam);
    } catch (const std::exception& exc) {
      MP_RAISE(std::string("Adding constraint of type '") +
                             Constraint::GetTypeName() + "' to " +
                             Backend::GetTypeName() + std::string(": ") +
                             exc.what());
    }
  }

  /// Log constraint group
  void LogConstraintGroup(
      BasicFlatModelAPI& be) override {
    auto cg = GetConstraintGroup(be);
    if (cg>=0 && GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["CON_GROUP"] = ConGroupName(cg);
        jw["CON_GROUP_index"] = cg;
      }
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
    }
  }


protected:
  /// Retrieve the Converter, const
  const Converter& GetConverter() const { return cvt_; }
  /// Retrieve the Converter
  Converter& GetConverter() { return cvt_; }

  /// Check constraint index
  bool check_index(int i) const { return i>=0 && i<(int)cons_.size(); }

  /// Container for a single constraint
  class Container {
  public:
    Container(int d, Constraint&& c) noexcept
      : con_(std::move(c)), depth_(d) { }

    /// Depth in redef tree
    int GetDepth() const { return depth_; }

    /// Bridged (reformulated or just unused.)
    /// If only reformulated, can still be checked
    /// for solution correctness.
    bool IsBridged() const { return is_bridged_; }
    /// Mark as bridged
    void MarkAsBridged() { is_bridged_=true; }

    /// Unused (should not be checked)
    bool IsUnused() const { return is_unused_; }
    /// Mark as unused
    void MarkAsUnused() {
      MarkAsBridged();
      is_unused_=true;
    }

    /// Get the flat constraint, const &
    const Constraint& GetCon() const { return con_.GetFlatConstraint(); }
    /// Get the flat constraint &
    Constraint& GetCon() { return con_.GetFlatConstraint(); }

  private:
    // Storing in the ExprWrapper,
    // so we can send (wrapper &) to ModelAPI::AddExpression().
    FlatExprType con_;
    int depth_ = 0;
    bool is_bridged_ = false;
    bool is_unused_ = false;
  };

	/// Convert all new constraints of this type
  bool ConvertAllFrom(int& i_last) {
    int i=i_last;
    const auto acceptanceLevel =
        GetChosenAcceptanceLevel();
    if (ConstraintAcceptanceLevel::NotAccepted == acceptanceLevel) {
      for ( ; ++i!=(int)cons_.size(); )
        if (!cons_[i].IsBridged())
          ConvertConstraint(cons_[i], i);
    }
    else if (ConstraintAcceptanceLevel::AcceptedButNotRecommended == acceptanceLevel) {
      for (; ++i != (int)cons_.size(); ) {
        if (!cons_[i].IsBridged()) {
          try {       // Try to convert all but allow failure
            ConvertConstraint(cons_[i], i);
          } catch (const ConstraintConversionGracefulFailure& ) {
            /// nothing
          } catch (const ConstraintConversionFailure& ccf) {
            GetConverter().AddWarning( ccf.key(), ccf.message() );
          }
        }
      }
    } else { // Recommended == acceptanceLevel &&
      for (; ++i != (int)cons_.size(); )
        if (!cons_[i].IsBridged() &&
            GetConverter().IfNeedsConversion(cons_[i].GetCon(), i))
          ConvertConstraint(cons_[i], i);
    }
    bool any_converted = i_last!=i-1;
    i_last = i-1;
    return any_converted;
  }

  void DoMarkForResultVars() {
    const auto eal        // expr only
        = GetChosenAcceptanceLevelEXPR();
    for (int i=0; i< (int)cons_.size(); ++i) {
      const auto& cnt = cons_[i];
      if (!cnt.IsBridged()) {      // Delegate actual logic to Converter
        const auto& con = cnt.GetCon();
        GetConverter().ConsiderMarkingResultVar(con, i, eal);
      }
    }
  }

  void DoMarkForArguments() {
    const auto eal        // expr only
        = GetChosenAcceptanceLevelEXPR();
    for (int i=0; i< (int)cons_.size(); ++i) {
      const auto& cnt = cons_[i];
      if (!cnt.IsBridged()) {      // Delegate actual logic to Converter
        const auto& con = cnt.GetCon();
        GetConverter().ConsiderMarkingArguments(con, i, eal);
      }
    }
  }

  void DoCvtWithExprs() { }

	/// Call Converter's RunConversion() and mark as "bridged".
  ///
	/// @param cnt the constraint container -
	/// actually redundant, as \a i is enough to find it. But for speed.
  /// @param i constraint index, needed for bridging
  void ConvertConstraint(Container& cnt, int i) {
    assert(!cnt.IsBridged());
    GetConverter().RunConversion(cnt.GetCon(), i, cnt.GetDepth());
    MarkAsBridged(cnt, i);
  }

  /// Mark item as reformulated
  void MarkAsBridged(Container& cnt, int ) {
		cnt.MarkAsBridged();
    ++n_bridged_or_unused_;
	}

  /// Mark item as unused
  void MarkAsUnused(Container& cnt, int ) {
    cnt.MarkAsUnused();
    ++n_bridged_or_unused_;
  }

protected:
  /// Export (last added) constraint
  void ExportConstraint(int i_con, const Container& cnt) {
    if (GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["index"] = i_con;
        if (*cnt.GetCon().name())
          jw["name"] = cnt.GetCon().name();
        jw["depth"] = cnt.GetDepth();
        WriteJSON(jw["data"], cnt.GetCon());
      }
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
    }
  }
  /// Export constraint status.
  /// This is called in the end,
  /// so printing the readable form.
  void ExportConStatus(int i_con, const Container& cnt,
                       const std::vector<std::string>* pvnam,
                       bool add2final) {
    if (GetLogger()) {
      fmt::MemoryWriter wrt;
      {
        MiniJSONWriter jw(wrt);
        jw["CON_TYPE"] = GetShortTypeName();
        jw["index"] = i_con;
        if (*cnt.GetCon().name()) {
          jw["name"] = cnt.GetCon().name();
          if (pvnam && pvnam->size()) {
            fmt::MemoryWriter pr;
            WriteFlatCon(pr, cnt.GetCon(), *pvnam);
            jw["printed"] = pr.c_str();
          }
        }
        jw["depth"] = cnt.GetDepth();
        jw["unused"] = (int)cnt.IsUnused();
        jw["bridged"] = (int)cnt.IsBridged();
        jw["final"] = (int)add2final;
      }
      wrt.write("\n");                     // EOL
      GetLogger()->Append(wrt);
    }
  }

public:
  /// Mark cons[\a i] as reformulated.
  /// Use index only.
  void MarkAsBridged(int i) override {
    MarkAsBridged(cons_.at(i), i);
	}

  /// Mark cons[\a i] as unused.
  /// Use index only.
  void MarkAsUnused(int i) override {
    MarkAsUnused(cons_.at(i), i);
  }

  /// Is constraint \a i unused?
  bool IsUnused(int i) const override {
    return cons_.at(i).IsUnused();
  }

  /// Copy names from ValueNodes
  void CopyNamesFromValueNodes() override {
    const auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      cons_[i].GetCon().SetName(vn[i].MakeCurrentName());
  }

  /// Copy names to ValueNodes
  void CopyNames2ValueNodes() {
    auto& vn = GetValueNode().GetStrVec();
    assert(vn.size()==cons_.size());
    for (auto i=vn.size(); i--; )
      vn[i] = std::string(cons_[i].GetCon().name());
  }

  /// ForEachActive().
  /// Deletes every constraint where fn() returns true.
	template <class Fn>
	void ForEachActive(Fn fn) {
		for (int i=0; i<(int)cons_.size(); ++i)
			if (!cons_[i].IsBridged())
        if (fn(cons_[i].GetCon(), i))
          MarkAsBridged(cons_[i], i);
	}

  /// Compute result for constraint \a i
  /// (for functional constraints).
  double
  ComputeValue(int i, const VarInfoRecomp& vir) override {
    assert(cons_[i].GetCon().GetResultVar() >= 0);
    return mp::ComputeValue(cons_[i].GetCon(), vir);
  }

  /// Compute violations for this constraint type.
  /// We do it for redefined (intermediate) ones too.
  void ComputeViolations(SolCheck& chk) override {
    if (cons_.size()) {
      auto& conviolmap =
          cons_.front().GetCon().IsLogical() ?
            chk.ConViolLog() :
            chk.ConViolAlg();
      const auto& x = chk.x_ext();
      ViolSummArray<3>* conviolarray {nullptr};
      for (int i=(int)cons_.size(); i--; ) {
        if (!cons_[i].IsUnused()) {
          int c_class = 0;    // class of this constraint
          if (!cons_[i].IsBridged())
            c_class |= 8;     // solver-side constraint
          if (!cons_[i].GetDepth())
            c_class |= 2;     // top-level
          if (!c_class)
            c_class = 4;      // intermediate
          if (c_class & chk.check_mode()) {
            auto viol = cons_[i].GetCon().ComputeViolation(x);
            auto cr = viol.Check(
                  chk.GetFeasTol(), chk.GetFeasTolRel());
            if (cr.first) {
              if (!conviolarray)
                conviolarray =         // lazy map access
                    &conviolmap[GetShortTypeName()];
              /// index==0,1,2: original, interm, solver-side
              /// If both orig and solver, report as orig
              int index = (c_class & 2) ? 0
                                        : (c_class & 8)
                                          ? 2 : 1;
              assert(index < (int)conviolarray->size());
              (*conviolarray)[index].CountViol(
                    viol, cr.second, cons_[i].GetCon().name());
            }
          }
        }
      }
    }
  }


protected:
  /// Add all non-converted items to ModelAPI.
  /// Export all constraints if desired.
  void AddAllUnbridged(BasicFlatModelAPI& be,
                       const std::vector<std::string>* pvnam) {
    int con_index=0;
    auto con_group = GetConstraintGroup(be);
    for (const auto& cont: cons_) {
      bool adding = !cont.IsBridged();
      if (adding) {
        static_cast<Backend&>(be).AddConstraint(cont.GetCon());
        GetConverter().GetCopyLink().
            AddEntry({
                       GetValueNode().Select(con_index),
                       GetConverter().GetValuePresolver().GetTargetNodes().
                         GetConValues()(con_group).Add()
                     });
      }
      ExportConStatus(con_index, cont, pvnam, adding);
      ++con_index;                      // increment index
    }
  }


private:
  Converter& cvt_;
  std::deque<Container> cons_;
  int i_cvt_last_ = -1;               // Last converted constraint.
  int n_bridged_or_unused_ = 0;       // Number of converted items,
                                      // they won't go to Backend
  const std::string desc_ {
    std::string("ConstraintKeeper< ") +
        Converter::GetTypeName() + ", " +
        Backend::GetTypeName() + ", " +
        Constraint::GetTypeName() + " >"};
};


////////////////////////////////////////////////////////////////////////////////////
/// Macros to define / access constraint keepers
/// Assume ConstraintManager as public parent

/// Use this to obtain a certain keeper, const
#define GET_CONST_CONSTRAINT_KEEPER(Constraint) \
  MPCD( GetConstraintKeeper((Constraint*)nullptr) )
/// Use this to obtain a certain keeper
#define GET_CONSTRAINT_KEEPER(Constraint) \
  MPD( GetConstraintKeeper((Constraint*)nullptr) )

/// Use this to obtain a certain map, const
#define GET_CONST_CONSTRAINT_MAP(Constraint) \
  MPCD( GetConstraintMap((Constraint*)nullptr) )
/// Use this to obtain a certain map
#define GET_CONSTRAINT_MAP(Constraint) \
  MPD( GetConstraintMap((Constraint*)nullptr) )

/// Define a constraint keeper
/// without a subexpression map.
/// @param optionNames: name (or a few, space-separated)
/// of the solver option(s) for acceptance of this
/// constraint
#define STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
private: \
  ConstraintKeeper<Impl, ModelAPI, Constraint> \
    CONSTRAINT_KEEPER_VAR(Constraint) \
      {*static_cast<Impl*>(this), \
       #Constraint, optionNames}; \
public: \
  const ConstraintKeeper<Impl, ModelAPI, Constraint>& \
  GetConstraintKeeper(const Constraint* ) const { \
    return CONSTRAINT_KEEPER_VAR(Constraint); \
  } \
  ConstraintKeeper<Impl, ModelAPI, Constraint>& \
  GetConstraintKeeper(Constraint* ) { \
    return CONSTRAINT_KEEPER_VAR(Constraint); \
  }

/// Define a constraint keeper
/// without a subexpression map.
/// Provide empty MapFind (returns -1) / MapInsert
#define STORE_CONSTRAINT_TYPE__NO_MAP( \
    Constraint, optionNames) \
  STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
  int MapFind__Impl(const Constraint& ) { return -1; } \
  bool MapInsert__Impl(const Constraint&, int ) \
    { return true; }


/// Define a constraint keeper
/// with a subexpression map.
/// The Converter storing the Constraint
/// should define MapFind / MapInsert accessing
/// the GET_(CONST_)CONSTRAINT_MAP(Constraint)
#define STORE_CONSTRAINT_TYPE__WITH_MAP( \
    Constraint, optionNames) \
  STORE_CONSTRAINT_TYPE__INTERNAL( \
    Constraint, optionNames) \
  STORE_CONSTRAINT_MAP(Constraint)

/// Internal use. Name of the constraint container
#define CONSTRAINT_KEEPER_VAR(Constraint) \
  ck__ ## Constraint ## _

/// Create constraint map. Normally internal use
#define STORE_CONSTRAINT_MAP(Constraint) \
  ConstraintMap<Constraint> CONSTRAINT_MAP_VAR(Constraint); \
  const ConstraintMap<Constraint>& \
  GetConstraintMap(Constraint* ) const { \
    return CONSTRAINT_MAP_VAR(Constraint); \
  } \
  ConstraintMap<Constraint>& \
  GetConstraintMap(Constraint* ) { \
    return CONSTRAINT_MAP_VAR(Constraint); \
  }

/// Internal use. Name of the constraint map
#define CONSTRAINT_MAP_VAR(Constraint) \
  map__ ## Constraint ## _


/////////////////////////////////////////////////////////////////////////////
/// Subexpression map
///
/// Indexes constraint location
template <class Constraint>
using ConstraintMap = std::unordered_map<
    std::reference_wrapper< const Constraint >, int >;

/// Subexpression maps for an expression require
/// operator==(refwrap<expr>, refwrap<expr>).
///
/// The below one is for CustomFunctionalConstraint<>
template <class Args, class Params, class NumOrLogic, class Id>
inline
bool operator==(std::reference_wrapper<
                  const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id> > c1,
                std::reference_wrapper<
                  const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id> > c2) {
  return c1.get().GetArguments() == c2.get().GetArguments() &&
      c1.get().GetParameters() == c2.get().GetParameters();
}

/// operator==(refwrap<ConditionalConstraint<> >)
template <class Con>
inline
bool operator==(std::reference_wrapper<
                  const ConditionalConstraint<Con> > c1,
                std::reference_wrapper<
                  const ConditionalConstraint<Con> > c2) {
  return c1.get().GetConstraint() == c2.get().GetConstraint();
}

/// operator==(LFC)
inline
bool operator==(std::reference_wrapper<
                  const LinearFunctionalConstraint > c1,
                std::reference_wrapper<
                  const LinearFunctionalConstraint > c2) {
  return c1.get().GetAffineExpr() == c2.get().GetAffineExpr();
}

/// operator==(QFC)
inline
bool operator==(std::reference_wrapper<
                  const QuadraticFunctionalConstraint > c1,
                std::reference_wrapper<
                  const QuadraticFunctionalConstraint > c2) {
  return c1.get().GetQuadExpr() == c2.get().GetQuadExpr();
}


////////////////////////////////////////////////////////////////////////////////////
/// Manage ConstraintKeepers for different constraint types
class ConstraintManager {
public:
  /// Add a new CKeeper with given conversion priority (smaller = sooner)
  void AddConstraintKeeper(BasicConstraintKeeper& ck, double priority) {
    con_keepers_.insert( { priority, ck } );
    ck.SetLogger(&*graph_exporter_app_);
  }

  /// This should be called after adding all constraint keepers
  void ConsiderAcceptanceOptions(
      BasicFlatConverter& cvt,
      const BasicFlatModelAPI& ma,
      Env& env) {
    for (auto& ck: con_keepers_)
      ck.second.ConsiderAcceptanceOptions(cvt, ma, env);
  }

  /// Convert all constraints (including any new appearing)
  void ConvertAllConstraints(BasicFlatConverter& cvt) {
    bool any_converted;
    do {
      any_converted = false;
      for (auto& ck: con_keepers_)
        any_converted = any_converted || ck.second.ConvertAllNewWith(cvt);
    } while (any_converted);
  }

  /// Mark which func cons can be expressions
  void MarkExprResultVars(BasicFlatConverter& cvt) {
    for (auto& ck: con_keepers_)
      ck.second.MarkExprResultVars(cvt);
  }

  /// Then, mark arguments of flat cons as vars
  void MarkArguments(BasicFlatConverter& cvt) {
    for (auto& ck: con_keepers_)
      ck.second.MarkArguments(cvt);
  }

  /// Convert to expression-based model
  void ConvertWithExpressions(BasicFlatConverter& cvt) {
    for (auto& ck: con_keepers_)
      ck.second.ConvertWithExpressions(cvt);
  }

  /// Fill counters of unbridged constraints
  void FillConstraintCounters(
      const BasicFlatModelAPI& mapi, FlatModelInfo& fmi) const {
    fmi.InitConstraintCount();
    for (const auto& ck: con_keepers_) {
      fmi.AddNumberOfConstraints(
        ck.second.GetTypeInfo(),
            ck.second.GetConstraintGroup(mapi),
            ck.second.GetNumberOfAddable());
    }
  }

  /// Copy names from ValueNodes
  void CopyNamesFromValueNodes() {
    for (const auto& ck: con_keepers_)
      ck.second.CopyNamesFromValueNodes();
  }

  /// Add all unbridged constraints to Backend
  void AddUnbridgedConstraintsToBackend(
      BasicFlatModelAPI& be,
      const std::vector<std::string>* pvnam=nullptr) const {
    for (const auto& ck: con_keepers_)
      ck.second.AddUnbridgedToBackend(be, pvnam);
  }

  /// Log constraint groups
  void LogConstraintGroups(
      BasicFlatModelAPI& be) const {
    for (const auto& ck: con_keepers_)
      ck.second.LogConstraintGroup(be);
  }

  /// Compute violations
  void ComputeViolations(SolCheck& chk) {
    for (const auto& ck: con_keepers_)
      ck.second.ComputeViolations(chk);
  }

  /// Retrieve file logger
  BasicFileAppender& GetFileAppender() const
  { return *graph_exporter_app_; }

private:
  std::multimap<double, BasicConstraintKeeper&> con_keepers_;
  /// Conversion graph exporter file appender
  std::unique_ptr<BasicFileAppender>
    graph_exporter_app_{MakeFileAppender()};
};

} // namespace mp

#endif // CONSTRAINT_KEEPER_H
