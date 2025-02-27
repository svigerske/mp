#ifndef MP2MIP_H
#define MP2MIP_H

#include <climits>

#include "mp/flat/converter.h"
#include "mp/flat/redef/MIP/redefs_mip_std.h"
#include "mp/flat/redef/encodings.h"

namespace mp {

/// MIPFlatConverter: converts flattened expressions for MIP
template <class Impl, class ModelAPI,
          class Model = FlatModel< > >
class MIPFlatConverter
    : public FlatConverter<Impl, ModelAPI, Model>
{
public:
  /// Class name for diagnostics
  static constexpr const char* GetTypeName() { return "MIPFlatConverter"; }

  /// BaseConverter typedef
  using BaseConverter = FlatConverter<Impl, ModelAPI, Model>;

  /// Constructor
  MIPFlatConverter(Env& e) : BaseConverter(e) {  }


  ///////////////////// MIP CONSTRAINT CONVERTERS ///////////////////////
  /// Reuse default (empty) ones
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter );

  /// Abs
  INSTALL_ITEM_CONVERTER(AbsConverter_MIP)
  /// And
  INSTALL_ITEM_CONVERTER(AndConverter_MIP)
  /// AllDiff
  INSTALL_ITEM_CONVERTER(AllDiffConverter_MIP)
  /// Complementarity linear
  INSTALL_ITEM_CONVERTER(ComplCvtLin_MIP)
  /// Complementarity quadratic
  INSTALL_ITEM_CONVERTER(ComplCvtQuad_MIP)
  /// Count
  INSTALL_ITEM_CONVERTER(CountConverter_MIP)
  /// Div
  INSTALL_ITEM_CONVERTER(DivConverter_MIP)
  /// CondLinEQ
  INSTALL_ITEM_CONVERTER(CondLinEQConverter_MIP)
  /// CondQuadEQ
  INSTALL_ITEM_CONVERTER(CondQuadEQConverter_MIP)
  /// CondLinLE
  INSTALL_ITEM_CONVERTER(CondLinLEConverter_MIP)
  /// CondQuadLE
  INSTALL_ITEM_CONVERTER(CondQuadLEConverter_MIP)
  /// CondLinLT
  INSTALL_ITEM_CONVERTER(CondLinLTConverter_MIP)
  /// CondQuadLT
  INSTALL_ITEM_CONVERTER(CondQuadLTConverter_MIP)
  /// CondLinGE
  INSTALL_ITEM_CONVERTER(CondLinGEConverter_MIP)
  /// CondQuadGE
  INSTALL_ITEM_CONVERTER(CondQuadGEConverter_MIP)
  /// CondLinGT
  INSTALL_ITEM_CONVERTER(CondLinGTConverter_MIP)
  /// CondQuadGT
  INSTALL_ITEM_CONVERTER(CondQuadGTConverter_MIP)
  /// IfThenElse
  INSTALL_ITEM_CONVERTER(IfThenElseConverter_MIP)
  /// Implication
  INSTALL_ITEM_CONVERTER(ImplicationConverter_MIP)
  /// ImplLE
  INSTALL_ITEM_CONVERTER(IndicatorLinLEConverter_MIP)
  /// ImplEQ
  INSTALL_ITEM_CONVERTER(IndicatorLinEQConverter_MIP)
  /// ImplGE
  INSTALL_ITEM_CONVERTER(IndicatorLinGEConverter_MIP)
  /// ImplQuadLE
  INSTALL_ITEM_CONVERTER(IndicatorQuadLEConverter_MIP)
  /// ImplQuadEQ
  INSTALL_ITEM_CONVERTER(IndicatorQuadEQConverter_MIP)
  /// ImplQuadGE
  INSTALL_ITEM_CONVERTER(IndicatorQuadGEConverter_MIP)
  /// Min
  INSTALL_ITEM_CONVERTER(MinConverter_MIP)
  /// Max
  INSTALL_ITEM_CONVERTER(MaxConverter_MIP)
  /// Not
  INSTALL_ITEM_CONVERTER(NotConverter_MIP)
  /// NumberofConst
  INSTALL_ITEM_CONVERTER(NumberofConstConverter_MIP)
  /// NumberofVar
  INSTALL_ITEM_CONVERTER(NumberofVarConverter_MIP)
  /// Or
  INSTALL_ITEM_CONVERTER(OrConverter_MIP)
  /// PLConstraint
  INSTALL_ITEM_CONVERTER(PLConverter_MIP)
  /// Pow
  INSTALL_ITEM_CONVERTER(PowConverter_MIP)
  /// PowConstExp
  INSTALL_ITEM_CONVERTER(PowConstExponentConverter_MIP)
  /// SOS2
  INSTALL_ITEM_CONVERTER(SOS2Converter_MIP)

  /// Smooth nonlinear

  /// Multiplication
  INSTALL_ITEM_CONVERTER(MulCvtLE_MIP)
  INSTALL_ITEM_CONVERTER(MulCvtEQ_MIP)
  INSTALL_ITEM_CONVERTER(MulCvtGE_MIP)

  /// Exp
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Exp)
  /// ExpA
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_ExpA)
  /// Log
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Log)
  /// LogA
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_LogA)
  /// Sin
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Sin)
  /// Cos
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Cos)
  /// Tan
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Tan)
  /// ASin
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Asin)
  /// ACos
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Acos)
  /// ATan
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Atan)
  /// Sinh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Sinh)
  /// Cosh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Cosh)
  /// Tanh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Tanh)
  /// ASinh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Asinh)
  /// ACosh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Acosh)
  /// ATanh
  INSTALL_ITEM_CONVERTER(FuncConConverter_MIP_Atanh)


  /// Strict comparison tolerance.
  /// Need a big eps to avoid misinterpretation,
  /// at least the solver's feasibility tolerance
  double ComparisonEps(int var) const {
    return MPCD(is_var_integer(var)) ? 1.0 : cmpEpsContinuous();
  }
  /// Strict comparison tolerance
  double ComparisonEps(var::Type vartype) const {
    return var::INTEGER==vartype ? 1.0 : cmpEpsContinuous();
  }


  ///////////////////////////////////////////////////////////////////////
  /////////////////////////// MIP CONSTRAINT MAPS ///////////////////////
public:
  USE_BASE_MAP_FINDERS( BaseConverter )

  /// Specialize MapFind for CondLinConEQ.
  /// We distinguish the case var==const
  int MapFind(const CondLinConEQ& eq0c) {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first)                    // only var==const comparisons
      return MapFind__VarConstCmp(isVCC.second.first, isVCC.second.second);
    return MPD( MapFind__Impl(eq0c) );
  }

  /// Specialize MapInsert for CondLiConEQ
  bool MapInsert(const CondLinConEQ& eq0c, int i) {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first)                    // only var==const comparisons
      return MapInsert__VarConstCmp(isVCC.second.first, isVCC.second.second, i);
    return MPD( MapInsert__Impl(eq0c, i) );
  }

  /// Convert MIP-specific maps
  void ConvertMaps() {
    BaseConverter::ConvertMaps();
    MP_DISPATCH( ConvertEqVarConstMaps() );
  }

  /// Whether might use equality encoding for \a var.
  /// @return false when UEnc apriori not applicable
  bool IfMightUseEqualityEncodingForVar(int var) const {
    if (n_map_cvt_cycles_)   // we've started map conversion
      return false;
    if (IsUEncAprioriExcluded(var))
      return false;
    if (!MPCD(is_var_integer(var))) {
      ExcludeUEncApriori(var);
      return false;
    }
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    bool fMight = (lb_dbl>std::numeric_limits<int>::min() &&
                   ub_dbl<std::numeric_limits<int>::max()) &&
                  (ub_dbl-lb_dbl <= 10000000);
    if (!fMight)
      ExcludeUEncApriori(var);
    return fMight;
  }

  /// Mark variable as "apriori no UEnc"
  void ExcludeUEncApriori(int var) const
  { set_vars_noUEnc_apriori_.insert(var); }
  /// Is var apriori not for UEnc?
  bool IsUEncAprioriExcluded(int var) const {
    return set_vars_noUEnc_apriori_.end() != set_vars_noUEnc_apriori_.find(var);
  }

  /// Obtain extended column \a k of ZZI encoding C^r
  std::vector<double>
  GetZZIExtendedColumn(int r, int k, int v0, int v1) {
    return GetExtendedColumn(*p_zzi_, r, k, v0, v1);
  }


private:
  /// For a single variable, map its equality comparisons
  /// for the comparison value (double), map the EQ0Constraint index
  using SingleVarEqConstMap = std::unordered_map<double, int>;
  /// A map keeping such maps for certain variables
  using VarsEqConstMap = std::unordered_map<int, SingleVarEqConstMap>;
  /// A set for variables not subject to UEncoding apriori
  using AprioriNoUEncSet = std::unordered_set<int>;

  VarsEqConstMap map_vars_eq_const_;
  mutable AprioriNoUEncSet set_vars_noUEnc_apriori_;
  int n_map_cvt_cycles_ = 0;      // number of cycles for map conversion

  P_ZZI_Encoding p_zzi_ { MakeZZIEncoding() };


protected:
  /// Return index of CondLinEQ for this \a val,
  /// if saved before
  int MapFind__VarConstCmp(int var, double val) {
    auto itVar = map_vars_eq_const_.find(var);
    if (map_vars_eq_const_.end() != itVar) {
      auto itCmp = itVar->second.find( val );
      if (itVar->second.end() != itCmp)
        /// Make sure we store the comparisons in CondConEQ's
        return itCmp->second;
    }
    return -1;
  }

  bool MapInsert__VarConstCmp(int var, double val, int i) {
    auto result = map_vars_eq_const_[var].
        insert( { val, i } );
    return result.second;
  }

  /// Result of IsVarConstCmp(): var, const
  using VarConstCmp = std::pair<int, double>;

  /// Check if \a con is a conditional var==const
  static std::pair<bool, VarConstCmp>
  IsVarConstCmp(const CondLinConEQ& con) {
    const auto& linEQ = con.GetConstraint();
    if (1==linEQ.size()) {
      assert(1.0==linEQ.coef(0));
      return { true, { linEQ.var(0), linEQ.rhs() } };
    }
    return { false, {} };
  }

  /// Convert the var-to-const comparisons left until now
  void ConvertEqVarConstMaps() {
    if (n_map_cvt_cycles_++)
      MP_RAISE("Repeated map conversion cycle");
    for (const auto& m: map_vars_eq_const_) {
      if (!IsUEncAprioriExcluded(m.first)) {  // else it's all done
        if (!WentWithoutEqEncForVar(m.first, m.second)) {
          ConvertEqVarConstMap(m.first, m.second);
          MarkUnusedConditionals(m.first, m.second);
        }
      }
    }
    // 2. Initiate possible conversions into big-M's
    MPD( ConvertAllConstraints() );
  }   // But why the map then?

  /// Aposteriori decide & reformulate cond equalities
  /// without unary encoding?
  bool WentWithoutEqEncForVar(
      int var, const SingleVarEqConstMap& map) {
    if (DontNeedEqEncForVar(var, map)) {
      GoWithoutEqEnc(var, map);
      return true;
    }
    return false;
  }

  /// @return true if don't need UEnc for var _aposteriori_
  bool DontNeedEqEncForVar(
      int var, const SingleVarEqConstMap& map) {
    if (IsUEncAprioriExcluded(var))
      return true;
    const auto& ck = GET_CONSTRAINT_KEEPER(CondLinConEQ);
    if (!ck.IfConsiderConversion())
      return true;
    double dom_rng = MPCD(ub(var))-MPCD(lb(var))+1;
    if (options_.NoUEncPosCtxRatio_ * map.size()   // Use UEnc if >
        > dom_rng)
      return false;
    int nNegCtx = 0;
    for (const auto& el: map) {
      const auto& con = ck.GetConstraint(el.second);
      if (con.GetContext().HasNegative())
        ++nNegCtx;
    }  // When up to 1 value in negative ctx, allow indicators.
    // Example. x is conditionally equated to a single value 5:
    // x==5 ==> ...
    return nNegCtx <= options_.NoUEncNegCtxMax_;
  }

  /// Manually convert all comparisons for this variable
  void GoWithoutEqEnc(int , const SingleVarEqConstMap& map) {
    auto& ck = GET_CONSTRAINT_KEEPER(CondLinConEQ);
    if (!ck.IfConsiderConversion())
      return;
    // 1. Convert the ConLinEq's into indicators.
    // Make sure IfMightUseEqualityEncoding() returns false.
    for (const auto& el: map) {
      ck.ConvertConstraint(el.second);
    }
    // 2. Initiate possible conversions into big-M's
    //    - done later for all.
  }

  void ConvertEqVarConstMap(int var, const SingleVarEqConstMap& map) {
    CreateUnaryEncoding(var, map);
  }

  /// Create unary encoding once we go for it.
  /// Adds a dummy UnaryEncodingConstraint for graph linking.
  /// Another design would be to have an "EqualityEncodingConstraint"
  /// and move all conversions there (indicators vs unary encoding.)
  /// That would partially automate linking.
  /// A benefit would be if some solver later supports them natively
  /// (IloDistribute?), now they are always converted below.
  void CreateUnaryEncoding(int var, const SingleVarEqConstMap& map) {
    // Create links from CondLinEQ's into the dummy UEnc
    auto valnode = MPD( AddConstraint( UnaryEncodingConstraint{{var}} ) );
    auto& ck = GET_CONSTRAINT_KEEPER(CondLinConEQ);
    for (const auto& veq : map) {
      MPD( GetMany2OneLink().AddEntry(
        { ck.SelectValueNodeRange(veq.second), valnode } ));
    }
    // Start AutoLinking from the dummy UEnc into its reformulation
    pre::AutoLinkScope<Impl> auto_link_scope{
      *(Impl*)this, valnode
    };
    // Create UEncoding
    const Model& model = MP_DISPATCH( GetModel() );
    if (!model.is_integer_var(var))
      MP_RAISE("Unary encoding: comparing non-integer variables not implemented");
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    if (lb_dbl==this->MinusInfty() || ub_dbl==this->Infty())
      MP_RAISE("Equality-comparing unbounded variables not implemented");
    if (lb_dbl<std::numeric_limits<int>::min() || ub_dbl>std::numeric_limits<int>::max())
      MP_RAISE("Equality-comparing variables with domain out of integer range not implemented");
    const int lb = (int)std::round(lb_dbl);
    const int ub = (int)std::round(ub_dbl);
    std::vector<int> unaryEncVars(ub-lb+1);
    int nTaken=0;
    for (int v=lb; v!=ub+1; ++v) {  // Running thru the map faster?
      auto itV = map.find(v);
      if (map.end() != itV) {
        ++nTaken;
        unaryEncVars[v-lb] =
          GET_CONSTRAINT_KEEPER(CondLinConEQ).GetResultVar(itV->second);
      } else {
        unaryEncVars[v-lb] = int( this->AddVar(0.0, 1.0, var::INTEGER) );
      }
    }
    assert(map.size()==(size_t)nTaken);
    std::vector<double> coefs(ub-lb+1, 1.0);
    this->AddConstraint(LinConEQ({coefs, unaryEncVars}, 1.0));
    unaryEncVars.push_back(var);
    for (int v=lb; v!=ub+1; ++v) {
      coefs[v-lb] = v;
    }
    coefs.push_back(-1.0);
    this->AddConstraint(LinConEQ({coefs, unaryEncVars}, 0.0));
  }

  /// Mark unused CondLinEQ's which were replaced by UEnc
  void MarkUnusedConditionals(int var, const SingleVarEqConstMap& map) {
    auto& ck = GET_CONSTRAINT_KEEPER(CondLinConEQ);
    for (const auto& el: map) {
      assert(!ck.IsRedundant(el.second));
      ck.MarkAsBridged(el.second);
    }
  }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////           ///
public:
  /// Init MIPFlatConverter options
  void InitOptions() {
    BaseConverter::InitOptions();
    InitOwnOptions();
  }

public:
  double cmpEpsContinuous() const { return options_.cmpEps_; }
  double bigMDefault() const { return options_.bigM_default_; }
  double PLApproxRelTol() const { return options_.PLApproxRelTol_; }
  double PLApproxDomain() const { return options_.PLApproxDomain_; }

private:
  struct Options {
    double cmpEps_ { 1e-4 };
    double bigM_default_ { -1 };
    double PLApproxRelTol_ { 1e-2 };
    double PLApproxDomain_ { 1e6 };
    double NoUEncPosCtxRatio_ { 0.0 };
    int NoUEncNegCtxMax_ { 1 };
  };
  Options options_;

  void InitOwnOptions() {
    this->GetEnv().AddOption("cvt:mip:eps cvt:cmp:eps cmp:eps",
                       "Tolerance for strict comparison of continuous variables for MIP. "
                             "Applies to <, >, and != operators. "
                       "Also applies to negation of conditional comparisons: "
                       "b==1 <==> x<=5 means that with b==0, x>=5+eps. "
                             "Default: 1e-4.",
                       options_.cmpEps_, 0.0, 1e100);
    this->GetEnv().AddOption("cvt:bigM cvt:bigm cvt:mip:bigM cvt:mip:bigm",
                       "Default value of big-M for linearization of logical constraints. "
                       "Not used by default. Use with care (prefer tight bounds). "
                       "Should be smaller than (1.0 / [integrality tolerance])",
                       options_.bigM_default_, -1.0, 1e100);
    this->GetEnv().AddOption("cvt:plapprox:reltol plapprox:reltol plapproxreltol",
                       "Relative tolerance for piecewise-linear approximation. Default 0.01.",
                       options_.PLApproxRelTol_, 0.0, 1e100);
    this->GetEnv().AddOption("cvt:plapprox:domain plapprox:domain plapproxdomain",
                       "For piecewise-linear approximated functions, both arguments and result "
                       "are bounded to +-[pladomain]. Default 1e6.",
                       options_.PLApproxDomain_, 0.0, 1e100);
    this->GetEnv().AddOption("cvt:uenc:ratio uenc:ratio",
                       "Min ratio (ub-lb)/Nvalues to skip unary encoding "
                       "for a variable x, where Nvalues is the number of constants "
                       "used in conditional comparisons x==const. Instead, "
                       "indicator constraints (or big-Ms) are used, if "
                       "uenc:negctx also applies. Default 0.",
                       options_.NoUEncPosCtxRatio_, 0.0, 1e100);
    this->GetEnv().AddOption("cvt:uenc:negctx:max uenc:negctx:max uenc:negctx",
                       "If cvt:uenc:ratio applies, max number of constants "
                       "in comparisons x==const in negative context "
                       "(equivalently, x!=const in positive context) to skip "
                       "UEnc(x). Default 1.",
                       options_.NoUEncNegCtxMax_, 0, INT_MAX);
  }
};

} // namespace mp

#endif // MP2MIP_H
