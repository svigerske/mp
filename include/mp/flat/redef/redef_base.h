#ifndef REDEF_BASE_H
#define REDEF_BASE_H

#include "mp/format.h"

namespace mp {

/// A base class for specific Item Converters, such as
/// individual constraint converters.
/// @param ModelConverter
template <class ModelConverter>
class BasicItemConverter {
public:
  /// Constructor
  BasicItemConverter(ModelConverter& mc) : mdl_cvt_(mc) { }

  /// Generic check whether the constraint
  /// needs to be converted despite being recommended for
  /// acceptance by ModelAPI.
  /// Example: PowConstExpConstraint(x, ...) with lb(x)<0.
  /// Default: false.
  template <class ItemType>
  bool IfNeedsConversion(const ItemType& , int ) {
    return false;
  }

  /// Generic check whether the constraint
  /// needs to be skipped from conversion,
  /// despite being not accepted by ModelAPI.
  /// Example: CondLinEQ(x == const).
  /// Default: false.
  template <class ItemType>
  bool IfDelayConversion(const ItemType& , int ) {
    return false;
  }

  /// Access const ModelConverter
  const ModelConverter& GetMC() const { return mdl_cvt_; }
  /// Access ModelConverter
  ModelConverter& GetMC() { return mdl_cvt_; }

private:
  ModelConverter& mdl_cvt_;
};


/// A functional constraint converter
/// @param Impl: the final converter
/// @param ModelConverter
template <class Impl, class ModelConverter>
class BasicFuncConstrCvt :
    public BasicItemConverter<ModelConverter> {
public:
  /// Base class
  using Base = BasicItemConverter<ModelConverter>;
  /// Constructor
  BasicFuncConstrCvt(ModelConverter& mc) : Base(mc) { }

  /// Generic Convert(), distinguishes context & result variable.
  /// Distinguish positive and negative context because
  /// if only one of these is necessary, significant model reduction
  /// can be achieved in certain cases.
  ///
  /// Responsible for adding presolve links, if any
  /// @param item: the item to be converted
  /// @param i: item index, used to create a presolve link
  ///
  /// The Impl can reimplement this
  template <class ItemType>
  void Convert(const ItemType& item, int i) {
    auto ctx = item.GetContext();
    assert(!ctx.IsNone());
    auto rv = item.GetResultVar();
    auto bnd00 = item.GetAprioriBounds();
    if ( ctx.HasNegative() &&
         GetMC().lb(rv) < bnd00.second ) {  // Need the negative direction
      MPD( ConvertCtxNeg(item, i) );
    }
    if ( ctx.HasPositive() &&
         GetMC().ub(rv) > bnd00.first ) {   // Need the positive direction
      MPD( ConvertCtxPos(item, i) );
    }
  }

  /// Convert in negative context
  template <class ItemType>
  void ConvertCtxNeg(const ItemType& item, int ) {
    MP_RAISE( fmt::format(
                "Conversion of '{}' in negative context "
                "not implemented", item.GetTypeName() ) );
  }
  /// Convert in positive context
  template <class ItemType>
  void ConvertCtxPos(const ItemType& item, int ) {
    MP_RAISE( fmt::format(
                "Conversion of '{}' in positive context "
                "not implemented", item.GetTypeName() ) );
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// To be used by descendants of BasicConverter
#define GET_CONSTRAINT_VALUE_NODE(con_type) \
  this->GetMC().GetValueNode((con_type*)nullptr)

/// In the ModelConverter: to use a specific item_cvt_type<>
/// Assumes Impl is the final ModelConverter type
#define INSTALL_ITEM_CONVERTER(item_cvt_type) \
  item_cvt_type<Impl> item_cvt__ ## item_cvt_type ## _ \
   { *static_cast<Impl*>(this) }; \
  bool IfHasCvt_impl(const typename \
      item_cvt_type<Impl>::ItemType*) { \
    return true; \
  } \
  bool IfNeedsCvt_impl(const typename \
      item_cvt_type<Impl>::ItemType& con, int i) { \
    return item_cvt__ ## item_cvt_type ## _ . \
      IfNeedsConversion(con, i); \
  } \
  bool IfDelayCvt_impl(const typename \
                     item_cvt_type<Impl>::ItemType& con, int i) { \
      return item_cvt__ ## item_cvt_type ## _ . \
      IfDelayConversion(con, i); \
  } \
  void Convert(const typename \
      item_cvt_type<Impl>::ItemType& con, int i) { \
    item_cvt__ ## item_cvt_type ## _ . Convert(con, i); \
  }


/// Conversion failure helper
class ConstraintConversionFailure {
  const char *key_, *msg_;
public:
  ConstraintConversionFailure(const char* key, const char* msg) noexcept :
    key_(key), msg_(msg) { }
  /// Failure type, used to display infos about failures
  const char* key() const { return key_; }
  /// Detailed message, should help improve model
  const char* message() const { return msg_; }
};


/// Graceful constraint conversion failure - no warnings
class ConstraintConversionGracefulFailure {
public:
  ConstraintConversionGracefulFailure() { }
};


} // namespace mp

#endif // REDEF_BASE_H
