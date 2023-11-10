/**
 * "Implementation" of NL Writer utils:
 * C++ classes "wrapping" them
 * in order to interface for NLWriter2, NLReader2
 */
#ifndef NLWRITER2MISCCIMPL_H
#define NLWRITER2MISCCIMPL_H

#include "api/c/nl-writer2-misc-c.h"

#include "mp/nl-writer2-misc.h"

namespace mp {

/// Wrap NLUtils_C into a C++ class,
/// in order to interface it for NLWriter2, NLReader2
class NLUtils_C_Impl
    : public NLUtils_C {
public:
  /// Construct
  NLUtils_C_Impl(NLUtils_C* pu)
    : p_nlu_c_(pu) { }

private:
  NLUtils_C* p_nlu_c_ {nullptr};
};

}  // namespace mp

#endif // NLWRITER2MISCCIMPL_H
