#ifndef MUL_H
#define MUL_H

#include "mp/error.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Convert all quadratic terms in the body of a QC.
/// This is the most universal way,
/// in case some converter adds QC.
/// An alternative would be, when flattening, to replace
/// every multiplication by
/// a special new constraint, and convert only these.
template <class ModelConverter, int sens>
class QCConverter_MIP :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  QCConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = QuadConRhs<sens>;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  bool IfNeedsConversion(const ItemType& , int ) {
		return !GetMC().IfPassQuadCon();
  }

  /// Skip conversion?
  bool IfDelayConversion(const ItemType& , int ) {
    return
        GetMC().IfPassQuadCon()
        || GetMC().IfWantNLOutput();     // Assume QuadExpr accepted
  }

  /// Conversion
	void Convert(const ItemType& qc, int ) {
		LinearizeQPTerms(qc);
  }


protected:
  using Base::GetMC;

  void LinearizeQPTerms(const ItemType& qc) {
    const auto& body = qc.GetBody();
    // Copy lin terms
    auto lin_terms = body.GetLinTerms();
    // Convert quadratic terms
    const auto& qp_terms = body.GetQPTerms();
    for (int i=0; i<(int)qp_terms.size(); ++i) {
      auto c = qp_terms.coef(i);
      auto x = qp_terms.var1(i);
      auto y = qp_terms.var2(i);
      lin_terms.add( LinearizeQPTerm(c, x, y) );
    }
    // Sort linear body. AS_ROOT propagates context
    lin_terms.sort_terms();
    GetMC().AddConstraint_AS_ROOT( LinConRhs< sens >{
                             lin_terms, qc.GetRhsOrRange() } );
  }

  LinTerms LinearizeQPTerm(double c, int x, int y) {
    if (GetMC().is_binary_var(x) ||
        GetMC().is_binary_var(y))
      return LinearizeProductWithBinaryVar(c, x, y);
    return LinearizeGeneralProduct(c, x, y);
  }

  LinTerms LinearizeProductWithBinaryVar(double c, int x, int y) {
    LinTerms lt;
    MP_ASSERT_ALWAYS(GetMC().is_binary_var(x) ||
                     GetMC().is_binary_var(y),
                     "Can only convert product with a binary variable");
    bool is_x_bin = GetMC().is_binary_var(x);
    auto i_bin = is_x_bin ? x : y;       // index of the (chosen) binary var
    auto i_other = is_x_bin ? y : x;                  // the other var index
    auto x_ifthen = GetMC().AssignResultVar2Args(
          IfThenConstraint {{      // no context as of now
                                   i_bin,
                                   i_other,
                                   int( GetMC().MakeFixedVar(0.0) )
                            }} );
    lt.add_term(c, x_ifthen);
    return lt;
  }

  /// c*x*y -> 0.5*c((x+y)^2-x^2-y^2)
  LinTerms LinearizeGeneralProduct(double c, int x, int y) {
    LinTerms lt;
    auto x_plus_y = GetMC().AssignResultVar2Args(
          LinearFunctionalConstraint{{  // = x+y
                                        { {1.0, 1.0}, {x, y} }, 0.0 }} );
    auto x_plus_y_pow2 = GetMC().AssignResultVar2Args(
          PowConstExpConstraint{  // no context as of now
                          VarArray1{x_plus_y}, {2.0} } );
    lt.add_term(0.5*c, x_plus_y_pow2);
    auto x_pow2 = GetMC().AssignResultVar2Args(
          PowConstExpConstraint{  // no context as of now
                          VarArray1{x}, {2.0} } );
    lt.add_term(-0.5*c, x_pow2);
    auto y_pow2 = GetMC().AssignResultVar2Args(
          PowConstExpConstraint{  // no context as of now
                          VarArray1{y}, {2.0} } );
    lt.add_term(-0.5*c, y_pow2);
    return lt;
  }
};


/// Typedef MulCvt<QuadConLE>
template <class MC>
using MulCvtLE_MIP = QCConverter_MIP<MC, -1>;

/// Typedef MulCvt<QuadConEQ>
template <class MC>
using MulCvtEQ_MIP = QCConverter_MIP<MC,  0>;

/// Typedef MulCvt<QuadConGE>
template <class MC>
using MulCvtGE_MIP = QCConverter_MIP<MC,  1>;


} // namespace mp

#endif // MUL_H
