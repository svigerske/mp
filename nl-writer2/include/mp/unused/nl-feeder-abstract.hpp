#ifndef NLFEEDERABSTRACT_HPP
#define NLFEEDERABSTRACT_HPP

#include "mp/nl-feeder-abstract.h"


namespace mp {

/// Wrap any reference
template <class Type>
class WrapReference {
public:
  /// Construct
  WrapReference(Type& ref) : ref_(ref) { }
  /// Get const ref
  const Type& Get() const { return ref_; }
  /// Get ref
  Type& Get() { return ref_; }
private:
  Type& ref_;
};


/// Wrap given vector writer
template <class VW>
class SparseVectorWriterWrp
    : public BasicSparseVectorWriter<typename VW::value_type>,
      public WrapReference<VW> {
  /// Construct
  SparseVectorWriterWrp(VW& vw) : WrapReference<VW>(vw) { }
  /// Number of elements left to write
  int NLeft() const override { return this->Get().NLeft(); }
  /// Write entry
  void Write(int index, typename VW::value_type value) override
  { this->Get().Write(index, value); }
};


/// Wrap vector writer factory
template <class VWF>
class SparseVectorWriterFactoryWrp
    : public BasicSparseVectorWriterFactory<typename VWF::value_type>,
      public WrapReference<VWF> {
public:
  using T = typename VWF::value_type;
  /// Construct
  SparseVectorWriterFactoryWrp(VWF& vwf) : WrapReference<VWF>(vwf) { }
  /// Make the SparseVectorWriter.
  /// @note std::size_t without include:
  ///   https://stackoverflow.com/questions/36594569/which-header-should-i-include-for-size-t.
  BasicSparseVectorWriter<T>&
  MakeVectorWriter(decltype(sizeof(int)) nnz) override {
    wrt_ = this->Get().MakeVectorWriter(nnz);
    vww_ = {wrt_};
    return vww_;
  }
private:
  typename VWF::writer_type wrt_;
  SparseVectorWriterWrp<typename VWF::writer_type> vww_;
};

template <typename Impl, typename ExprType>
template <class ObjGradWriterFactory>
void NLFeederAbstract<Impl, ExprType>::
    FeedObjGradient(int i, ObjGradWriterFactory& ) {

}


/// STOPPED HERE


/// Dispatcher method
template <class ObjExprWriter>
void FeedObjExpression(int , ObjExprWriter& ew);


/// Dispatcher method
template <class DefVarWriterFactory>
void FeedDefinedVariables(int i, DefVarWriterFactory& );


/// Dispatcher method
template <class VarBoundsWriter>
void FeedVarBounds(VarBoundsWriter& );


/// Dispatcher method
template <class ConBoundsWriter>
void FeedConBounds(ConBoundsWriter& );


/// Dispatcher method
template <class ConLinearExprWriterFactory>
void FeedLinearConExpr(int i, ConLinearExprWriterFactory& );

/// Dispatcher method
template <class ConExprWriter>
void FeedConExpression(int , ConExprWriter& ew);


/// Dispatcher method
template <class ExprWriter>
void FeedExpr(Expr e, ExprWriter& );


/// Dispatcher method
template <class PLSOSWriter>
void FeedPLSOS(PLSOSWriter& );


/// Dispatcher method
template <class ColSizeWriter>
void FeedColumnSizes(ColSizeWriter& );


/// Dispatcher method
template <class IGWriter>
void FeedInitialGuesses(IGWriter& );

/// Dispatcher method
template <class IDGWriter>
void FeedInitialDualGuesses(IDGWriter& );


/// Dispatcher method
template <class SuffixWriterFactory>
void FeedSuffixes(SuffixWriterFactory& );


/// Dispatcher method
template <class RowObjNameWriter>
void FeedRowAndObjNames(RowObjNameWriter& wrt);

/// Dispatcher method
template <class DelRowNameWriter>
void FeedDelRowNames(DelRowNameWriter& );

/// Dispatcher method
template <class ColNameWriter>
void FeedColNames(ColNameWriter& );

/// Dispatcher method
template <class UnusedVarNameWriter>
void FeedUnusedVarNames(UnusedVarNameWriter& );

/// Dispatcher method
template <class FixedVarNameWriter>
void FeedFixedVarNames(FixedVarNameWriter& );

/// Dispatcher method
template <class ObjOffsetWriter>
void FeedObjAdj(ObjOffsetWriter& );


}  // namespace mp

#endif // NLFEEDERABSTRACT_HPP
