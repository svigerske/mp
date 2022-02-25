#ifndef MP2MIP_H
#define MP2MIP_H

#include "mp/flat/converter.h"

#include "mp/flat/redef/MIP/redefs_mip_std.h"

namespace mp {

/// MIPFlatConverter: converts flattened expressions for MIP
template <class Impl, class Backend,
          class Model = BasicFlatModel< > >
class MIPFlatConverter
    : public FlatConverter<Impl, Backend, Model>
{
public:
  static constexpr const char* name() { return "MIPFlatConverter"; };

public:
  using BaseConverter = FlatConverter<Impl, Backend, Model>;

public:
  static const char* GetConverterName() { return "MIPFlatConverter"; }
  MIPFlatConverter(Env& e) : BaseConverter(e) {  }

  ///////////////////// SPECIALIZED CONSTRAINT CONVERTERS //////////////////
  USE_BASE_CONSTRAINT_CONVERTERS( BaseConverter );        ///< reuse default ones

  double ComparisonEps(int var) const {
    return MPCD(is_var_integer(var)) ? 1.0 : 1e-6; // TODO param
  }
  double ComparisonEps(var::Type vartype) const {
    return var::INTEGER==vartype ? 1.0 : 1e-6; // TODO param
  }


  void Convert(const IndicatorConstraintLinLE& indc) {
    auto binvar=indc.get_binary_var();
    auto ae = indc.to_lhs_expr();
    auto bnds = MPD( ComputeBoundsAndType(ae) );
    ConvertImplicationLE(binvar, indc.get_binary_value(), bnds, std::move(ae));
  }

  /// (b==val ==> ae<=0)
  void ConvertImplicationLE(int b, int val,
                   const PreprocessInfoStd& bnds, AffExp ae) {
    /// TODO fail if lb>0 +report .iis if requested
    /// TODO skip if ub<0
    if (bnds.ub() >= this->PracticallyInfty())
      throw ConstraintConversionFailure( "IndicatorInfBound",
          "The redefinition of a (possibly auxiliary) indicator constraint failed"
          " so it had to be passed to the solver."
          " Provide tight bounds on variables entering logical expressions, "
          "or set acc:ind_le=2");
    if (val)            // left condition is b==1
      ae += {{bnds.ub(), b}, -bnds.ub()};
    else
      ae += {{-bnds.ub(), b}, 0.0};
    MP_DISPATCH( AddConstraint(LinConLE(     /// Big-M constraint
        (LinTerms&&)ae, -ae.constant_term() )) );
  }

  /// b==val ==> c'x==d
  void Convert(const IndicatorConstraintLinEQ& indc) {
    auto binvar=indc.get_binary_var();
    auto ae = indc.to_lhs_expr();
    auto bnds = MPD( ComputeBoundsAndType(ae) );
    ConvertImplicationLE(binvar, indc.get_binary_value(), bnds, ae);
    ae.negate();
    bnds.NegateBounds();
    ConvertImplicationLE(binvar, indc.get_binary_value(), bnds, std::move(ae));
  }

  void Convert(const ConjunctionConstraint& conj) {
    assert(!conj.GetContext().IsNone());
    if (conj.GetContext().HasPositive())
      ConvertImplied(conj);
    if (conj.GetContext().HasNegative())
      ConvertReverseImplied(conj);
  }

  void ConvertReverseImplied(const ConjunctionConstraint& conj) {
    const auto& args = conj.GetArguments();
    auto flags = args;
    flags.push_back(conj.GetResultVar());
    std::vector<double> ones(args.size(), 1.0); // res+n-1 >= sum(args) in CTX-
    ones.push_back(-1.0);
    MP_DISPATCH( AddConstraint(
                   LinConLE({ones, flags},
                               {(double)args.size()-1} )) );
  }

  void ConvertImplied(const ConjunctionConstraint& conj) {
    std::array<double, 2> coefs{-1.0, 1.0};
    std::array<int, 2> vars{-1, conj.GetResultVar()};
    for (auto arg: conj.GetArguments()) {       // res <= arg[i] in CTX+
      vars[0] = arg;
      MP_DISPATCH( AddConstraint(
                     LinConLE({coefs, vars}, {0.0} )) );
    }
  }

  void Convert(const DisjunctionConstraint& disj) {
    assert(!disj.GetContext().IsNone());
    if (disj.GetContext().HasPositive())
      ConvertImplied(disj);
    if (disj.GetContext().HasNegative())
      ConvertReverseImplied(disj);
  }

  void ConvertImplied(const DisjunctionConstraint& disj) {
    const auto& args = disj.GetArguments();
    auto flags = args;
    flags.push_back(disj.GetResultVar());
    std::vector<double> ones(args.size(), 1.0);  // res <= sum(args) in CTX+
    ones.push_back(-1.0);
    MP_DISPATCH( AddConstraint(
                   LinConGE({ones, flags}, {0.0} )) );
  }

  void ConvertReverseImplied(const DisjunctionConstraint& disj) {
    std::array<double, 2> coefs{1.0, -1.0};
    std::array<int, 2> vars{-1, disj.GetResultVar()};
    for (auto arg: disj.GetArguments()) {        // res >= arg[i] in CTX-
      vars[0] = arg;
      MP_DISPATCH( AddConstraint(
                     LinConLE({coefs, vars}, {0.0} )) );
    }
  }

  void Convert(const IfThenConstraint& itc) {
    assert(!itc.GetContext().IsNone());
    const auto& args = itc.GetArguments();
    if (!this->is_fixed(args[1]) || !this->is_fixed(args[2]))
      throw std::logic_error("MP2MIP: IfThen with variable then/else arguments not implemented");
    else
      ConvertIfThen_constantThenElse(itc);
  }

  void ConvertIfThen_constantThenElse(const IfThenConstraint& itc) {
    const auto& args = itc.GetArguments();
    assert((this->is_fixed(args[1]) && this->is_fixed(args[2])));
    const double const1 = this->fixed_value(args[1]);
    const double const2 = this->fixed_value(args[2]);
    this->AddConstraint( LinearFunctionalConstraint(
                           itc.GetResultVar(),
    { {{const1-const2}, {args[0]}}, const2 } ) );
  }

  //////////////////// NUMBEROF CONST ///////////////////////
  void Convert(const NumberofConstConstraint& nocc) {
    const auto& args = nocc.GetArguments();
    const double k = nocc.GetParameters()[0];
    std::vector<double> coefs(args.size()+1, 1.0);
    std::vector<int> flags(args.size()+1, nocc.GetResultVar());
    for (size_t ivar = 0; ivar < args.size(); ++ivar) {
      flags[ivar] = this->AssignResultVar2Args(   // flag = (args[i]==k)
            EQ0Constraint( { {{1.0}, {args[ivar]}}, -k } ) );
    }
    coefs.back() = -1.0;
    this->AddConstraint( LinConEQ( {coefs, flags}, {0.0} ) );
  }

  //////////////////// NUMBEROF VAR ///////////////////////
  /// Very basic, could be improved
  void Convert(const NumberofVarConstraint& novc) {
    const auto& args = novc.GetArguments();
    std::vector<double> coefs(args.size(), 1.0);
    coefs.front() = -1.0;
    std::vector<int> flags(args.size(), novc.GetResultVar());
    for (size_t ivar = 1; ivar < args.size(); ++ivar) {
      flags[ivar] = this->AssignResultVar2Args(   // flag = (args[i]==args[0])
            EQ0Constraint(
                 { { {1.0, -1.0}, {args[ivar], args[0]} }, 0.0 } ) );
    }
    this->AddConstraint( LinConEQ( {coefs, flags}, {0.0} ) );
  }

  void Convert(const CountConstraint& cc) {
    const auto& args = cc.GetArguments();
    std::vector<double> coefs(args.size()+1, 1.0);
    coefs.back() = -1.0;
    std::vector<int> flags(args.size()+1, cc.GetResultVar());
    for (size_t ivar = 0; ivar < args.size(); ++ivar) {
      flags[ivar] = args[ivar];
      /// Force booleanize
      /// Either reify !=0 or constrain to 0..1? TODO param? TODO warning?
      if (!MPD( is_binary_var(args[ivar]) )) {
        auto feq0 = this->AssignResultVar2Args(   // feq0 = (args[i]==0)
            EQ0Constraint( { {{1.0}, {args[ivar]}}, 0.0 } ) );
        flags[ivar] = this->AssignResultVar2Args(   // flag = (args[i]!=0)
            NotConstraint( {feq0} ));
      }
    }
    this->AddConstraint( LinConEQ( {coefs, flags}, 0.0 ) );
  }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////////// MAPS //////////////////////////
  ///
private:
  /// For a single variable, map its equality comparisons
  /// for the comparison value (double), map the EQ0Constraint index
  using SingleVarEqConstMap = std::unordered_map<double, int>;
  /// A map keeping such maps for certain variables
  using VarsEqConstMap = std::unordered_map<int, SingleVarEqConstMap>;

  VarsEqConstMap map_vars_eq_const_;

public:
  //////////////////////////////// CONSTRAINT MAPS ///////////////////////////////////////
  ///
  USE_BASE_MAP_FINDERS( BaseConverter )

  AbstractConstraintLocation MapFind(const EQ0Constraint& eq0c) {
    const auto isVCC = IsVarConstCmp( eq0c );
    if (isVCC.first) {                    // only var==const comparisons
      return MapFind__VarConstCmp(isVCC.second.first, isVCC.second.second);
    }
    return { };
  }

  AbstractConstraintLocation MapFind__VarConstCmp(int var, double val) {
    auto itVar = map_vars_eq_const_.find(var);
    if (map_vars_eq_const_.end() != itVar) {
      auto itCmp = itVar->second.find( val );
      if (itVar->second.end() != itCmp)
        /// Make sure we store the comparisons in EQ0Con's
        return {&GET_CONSTRAINT_KEEPER(EQ0Constraint), itCmp->second};
    }
    return { };
  }

  bool MapInsert(ConstraintLocation<EQ0Constraint> ck) {
    const auto isVCC = IsVarConstCmp( ck.GetConstraint() );
    if (isVCC.first) {                    // only var==const comparisons
      return MapInsert__VarConstCmp(isVCC.second.first, isVCC.second.second, ck);
    }
    return true;
  }

  bool MapInsert__VarConstCmp(int var, double val,
                              ConstraintLocation<EQ0Constraint> ck) {
    auto result = map_vars_eq_const_[var].
        insert( std::make_pair( val, ck.GetIndex() ) );
    return result.second;
  }

  /// Result of IsVarConstCmp(): var, const
  using VarConstCmp = std::pair<int, double>;

  /// Check if the eq0c is a var==const
  static std::pair<bool, VarConstCmp> IsVarConstCmp(const EQ0Constraint& con) {
    const AffExp& args = con.GetArguments();
    if (1==args.size()) {
      assert(1.0==args.coef(0));
      return { true, { args.var(0), -args.constant_term() } };
    }
    return { false, {} };
  }

  void ConvertMaps() {
    MP_DISPATCH( ConvertEqVarConstMaps() );
  }

  void ConvertEqVarConstMaps() {
    for (const auto& m: map_vars_eq_const_) {
      if (IfUseEqualityEncodingForVar(m.first))
        ConvertEqVarConstMap(m.first, m.second);
    } // Otherwise, indicators / big-Ms should have been applied
  }

  bool IfUseEqualityEncodingForVar(int var) const {
    if (!MPCD(is_var_integer(var)))
      return false;
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    return (lb_dbl>std::numeric_limits<int>::min() && // TODO enable outside
            ub_dbl<std::numeric_limits<int>::max()) &&
        (ub_dbl-lb_dbl <= 100000);      // TODO param
  }

  void ConvertEqVarConstMap(int var, const SingleVarEqConstMap& map) {
    CreateUnaryEncoding(var, map);
  }

  void CreateUnaryEncoding(int var,  const SingleVarEqConstMap& map) {
    const Model& model = MP_DISPATCH( GetModel() );
    if (!model.is_integer_var(var))
      throw std::logic_error("MP2MIP: Equality-comparing non-integer variables not implemented");
    const auto lb_dbl = this->lb(var);
    const auto ub_dbl = this->ub(var);
    if (lb_dbl==this->MinusInfty() || ub_dbl==this->Infty())
      throw std::logic_error("MP2MIP: Equality-comparing unbounded variables not implemented");
    if (lb_dbl<std::numeric_limits<int>::min() || ub_dbl>std::numeric_limits<int>::max())
      throw std::logic_error("MP2MIP: Equality-comparing variables with domain out of integer range not implemented");
    const int lb = (int)std::round(lb_dbl);
    const int ub = (int)std::round(ub_dbl);
    std::vector<int> unaryEncVars(ub-lb+1);
    int nTaken=0;
    for (int v=lb; v!=ub+1; ++v) { // TODO run thru the map first
      auto itV = map.find(v);
      if (map.end() != itV) {
        ++nTaken;
        unaryEncVars[v-lb] =
          GET_CONSTRAINT_KEEPER(EQ0Constraint).GetResultVar(itV->second);
      } else {
        unaryEncVars[v-lb] = this->AddVar(0.0, 1.0, var::INTEGER);
      }
    }
    assert(map.size()==(size_t)nTaken);
    std::vector<double> coefs(ub_dbl-lb_dbl+1, 1.0);
    this->AddConstraint(LinConEQ({coefs, unaryEncVars}, 1.0));
    unaryEncVars.push_back(var);
    for (int v=lb; v!=ub+1; ++v) {
      coefs[v-lb] = v;
    }
    coefs.push_back(-1.0);
    this->AddConstraint(LinConEQ({coefs, unaryEncVars}, 0.0));
  }

public:
  /////////////////////// CONSTRAINT CONVERTERS /////////////////////////

  /// Abs
  INSTALL_ITEM_CONVERTER(AbsConverter_MIP)
  /// AllDiff
  INSTALL_ITEM_CONVERTER(AllDiffConverter_MIP)
  /// EQ0
  INSTALL_ITEM_CONVERTER(EQ0Converter_MIP)
  /// LE0
  INSTALL_ITEM_CONVERTER(LE0Converter_MIP)
  /// Min
  INSTALL_ITEM_CONVERTER(MinConverter_MIP)
  /// Max
  INSTALL_ITEM_CONVERTER(MaxConverter_MIP)
  /// Not
  INSTALL_ITEM_CONVERTER(NotConverter_MIP)


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// OPTIONS /////////////////////////
  ///
public:
  void InitOptions() {
    BaseConverter::InitOptions();
  }

};

} // namespace mp

#endif // MP2MIP_H
