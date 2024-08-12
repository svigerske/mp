/*
 Basic expression-based model API definitions.

 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov <gleb@ampl.com>
*/
#ifndef MODEL_API_BASE_H
#define MODEL_API_BASE_H

#include "mp/flat/model_api_base.h"
#include "mp/flat/nl_expr/constr_nl.h"

namespace mp {

/// ModelAPIs handling expression trees should derive from
template <class ExprType=void*>
class BasicExprModelAPI
    :public BasicFlatModelAPI {
public:
  using Expr = ExprType;
  /// Placeholder for GetTypeName()
  static const char* GetTypeName()    { return "BasicExprModelAPI"; }

/// A ModelAPI accepting NL trees should declare this.
///
/// - NotAccepted: not compiled
/// - AcceptedButNotRecommended: compiled but off by default (option acc:_expr)
/// - Recommended: on by default
#define ACCEPT_EXPRESSION_INTERFACE(val) \
  static constexpr ExpressionAcceptanceLevel \
  ExpressionInterfaceAcceptanceLevel() { return ExpressionAcceptanceLevel::val; }

/// Reuse inherited names
  USE_BASE_CONSTRAINT_HANDLERS(BasicFlatModelAPI)

};

}  // namespace mp

#endif // MODEL_API_BASE_H
