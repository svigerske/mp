/*
 Optimization problem

 Copyright (C) 2014 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#include <map>

#include "mp/problem.h"
#include "mp/problem-builder.h"
#include "expr-writer.h"

namespace mp {

void LinearExpr::SortTerms() {
  std::map<int, double> var_coef_map;
  for (int i=0; i<num_terms(); ++i)
    if (0.0!=std::fabs(coef(i)))
      var_coef_map[var_index(i)] += coef(i);
  terms_.clear();
  for (const auto& vc: var_coef_map) {
    if (0.0!=std::fabs(vc.second))
      AddTerm(vc.first, vc.second);
  }
}

template <typename Alloc>
int BasicProblem<Alloc>::GetSuffixSize(suf::Kind kind) {
  std::size_t size = 0;
  switch (kind & suf::KIND_MASK) {
  default:
    MP_ASSERT(false, "invalid suffix kind");
    // Fall through.
  case suf::VAR:
    size = vars_.capacity();
    break;
  case suf::CON:
    size = algebraic_cons_.capacity() +
        logical_cons_.capacity();
    break;
  case suf::OBJ:
    size = linear_objs_.capacity();
    break;
  case suf::PROBLEM:
    size = 1;
    break;
  }
  return static_cast<int>(size);
}

template <typename Alloc>
typename BasicProblem<Alloc>::LinearObjBuilder BasicProblem<Alloc>::AddObj(
    obj::Type type, NumericExpr expr, int num_linear_terms) {
  MP_ASSERT(linear_objs_.size() < MP_MAX_PROBLEM_ITEMS, "too many objectives");
  is_obj_max_.push_back(type != obj::MIN);
  linear_objs_.push_back(LinearExpr());
  LinearExpr &linear_expr = linear_objs_.back();
  linear_expr.Reserve(num_linear_terms);
  if (expr)
    SetNonlinearObjExpr(static_cast<int>(linear_objs_.size() - 1), expr);
  return LinearObjBuilder(&linear_expr);
}

template <typename Alloc>
void BasicProblem<Alloc>::SetComplementarity(
    int con_index, int var_index, ComplInfo info) {
  MP_ASSERT(0 <= con_index && con_index < num_algebraic_cons(),
            "invalid index");
  MP_ASSERT(0 <= var_index && var_index < num_vars(), "invalid index");
  if (compl_vars_.size() <= static_cast<std::size_t>(con_index)) {
    compl_vars_.reserve(algebraic_cons_.capacity());
    compl_vars_.resize(algebraic_cons_.size());
  }
  compl_vars_[con_index] = var_index + 1u;
  AlgebraicConInfo &con = algebraic_cons_[con_index];
  con.lb = info.con_lb();
  con.ub = info.con_ub();
}

template <typename Alloc>
void BasicProblem<Alloc>::SetInfo(const NLProblemInfo &info) {
  vars_.reserve(info.num_vars);
  is_var_int_.reserve(info.num_vars);
  is_obj_max_.reserve(info.num_objs);
  linear_objs_.reserve(info.num_objs);
  if (info.num_nl_objs != 0)
    nonlinear_objs_.reserve(info.num_objs);
  algebraic_cons_.reserve(info.num_algebraic_cons);
  if (info.num_compl_conds != 0)
    compl_vars_.reserve(info.num_algebraic_cons);
  if (info.num_nl_cons != 0)
    nonlinear_cons_.reserve(info.num_algebraic_cons);
  logical_cons_.reserve(info.num_logical_cons);
  int num_common_exprs = info.num_common_exprs();
  linear_exprs_.reserve(num_common_exprs);
  nonlinear_exprs_.reserve(num_common_exprs);
  common_expr_positions_.reserve(num_common_exprs);
}

template <typename Params>
const std::string& BasicProblem<Params>::item_name(
    int i, std::vector<std::string>& names,
    int n, const char* stub,  int n2, const char* stub2) {
  assert(0<=i && i<n);
  if ((int)names.size()<n) {
    int sz0 = (int)names.size();
    names.resize(n);
    int k=sz0;
    auto gen_names = [&names](
        int k, int n, const char* stub, int k_sub=0) {
      for ( ; k<n; k++) {
        assert(stub);
        int l = std::strlen(stub);
        assert(l);
        bool fbr = '['==stub[l-1];
        names[k] = stub;
        names[k] += std::to_string(k-k_sub+1);
        names[k] += (fbr ? ']' : '_');
      }
    };
    if (n2>=0 && n2<n) {
      if (k<n2) {
        gen_names(k, n2, stub);
        gen_names(n2, n, stub2, n2);
      } else
        gen_names(k, n, stub2, n2);
    } else
      gen_names(k, n, stub);
  }
  assert(names[i].size());
  return names[i];
}

// Instantiations

template class BasicProblem< >;

template void ReadNLFile(fmt::CStringRef filename, Problem &p, int flags);

template
void ReadNLString(NLStringRef str, Problem &p, fmt::CStringRef name, int flags);

template
void WriteExpr<typename Problem::ExprTypes>(
fmt::Writer &w,
const LinearExpr &linear, NumericExpr nonlinear,
Problem::VarNamer);

template
void WriteExpr<typename Problem::ExprTypes>
(fmt::Writer &w, LogicalExpr expr, Problem::VarNamer);

/// Write algebraic constraint.
template
void WriteAlgCon<typename Problem::ExprTypes>(
fmt::Writer &w, const typename Problem::MutAlgebraicCon &con,
Problem::VarNamer);

}  // namespace mp
