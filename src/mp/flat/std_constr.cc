#include <map>
#include <cfloat>
#include <cassert>

#include "mp/format.h"
#include "mp/common.h"

#include "mp/flat/expr_quadratic.h"
#include "mp/flat/obj_std.h"
#include "mp/flat/constr_keeper.h"
#include "mp/flat/model_info.hpp"

namespace mp {

/// @todo Keep consistent with the \a ConstraintGroups enum.
static const char* const congroup_names[]
= {
 "Default",
 "All",
 "Algebraic",
 "Linear",
 "Quadratic",
 "Conic",
 "General",
        "Nonlinear",
 "Piecewiselinear",
 "SOS",
 "SOS1",
 "SOS2",
 "Logical"
};

const char* ConGroupName(int cg) {
  assert(0<=cg
         && cg<int(sizeof (congroup_names)/sizeof(congroup_names[0])));
  return congroup_names[cg];
}

void LinTerms::sort_terms(bool force_sort) {
  if (1<size()
      || (1==size() && !coef(0))) {
    std::map<int, double> var_coef_map;
    for (size_t i=0; i<size(); ++i)
      if (0.0!=std::fabs(coefs_[i]))
        var_coef_map[vars_[i]] += coefs_[i];
    if (force_sort ||                    // force sorting for tests
        var_coef_map.size() < size()) {
      coefs_.clear();
      vars_.clear();
      for (const auto& vc: var_coef_map) {
        if (0.0!=std::fabs(vc.second)) {         // Need tolerance?
          coefs_.push_back(vc.second);
          vars_.push_back(vc.first);
        }
      }
    }
  }
}


void QuadTerms::sort_terms()  {
  if (1<size()
      || (1==size() && !coef(0))) {
    auto sort_pair = [](int a, int b) {
      return a<b ? std::pair<int, int>(a, b) : std::pair<int, int>(b, a);
    };
    std::map<std::pair<int, int>, double> var_coef_map;
    for (int i=0; i<size(); ++i)
      if (0.0!=std::fabs(coefs_[i]))
        var_coef_map[sort_pair(vars1_[i], vars2_[i])] += coefs_[i];
    if (true) {
      coefs_.clear();
      vars1_.clear();
      vars2_.clear();
      for (const auto& vc: var_coef_map) {
        if (0.0!=std::fabs(vc.second))         // Need tolerance?
          add_term(vc.second, vc.first.first, vc.first.second);
      }
    }
  }
}

const char*
BasicConstraintKeeper::GetShortTypeName() const {
  if (type_name_short_.empty()) {
    std::string acc_opt = GetAcceptanceOptionNames();
    assert(acc_opt.size());
    auto word_end = std::min(acc_opt.find(' '),
                             acc_opt.size());
    auto colon_pos = acc_opt.find(':');
    if (colon_pos>word_end)
      colon_pos = 0;
    type_name_short_ = acc_opt.substr(
        colon_pos, word_end-colon_pos);
    for (auto& c: type_name_short_)
      if (':'==c)
        c = '_';                // Markdown
    assert(type_name_short_.size());
  }
  return type_name_short_.c_str();
}


/// acceptance when constraint only
static const mp::OptionValueInfo values_con_acceptance[] = {
    { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
    { "1", "Accepted but automatic redefinition will be used where possible", 1},
    { "2", "Accepted natively and preferred", 2}
};

/// acceptance when expression only
static const mp::OptionValueInfo values_expr_acceptance[] = {
    { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
    { "3", "Accepted but automatic redefinition will be used where possible", 3},
    { "4", "Accepted natively and preferred", 4}
};

/// acceptance when both constraint and expression are possible
/// (ANY SOLVER DOING THIS? Ilog CP?)
static const mp::OptionValueInfo values_universal_acceptance[] = {
    { "0", "Not accepted natively, automatic redefinition will be attempted", 0},
    { "1", "Accepted as constraint but automatic redefinition will be used where possible", 1},
    { "2", "Accepted as constraint natively and preferred", 2},
    { "3", "Accepted as expression but automatic redefinition will be used where possible", 3},
    { "4", "Accepted as expression natively and preferred", 4}
};



void BasicConstraintKeeper::DoAddAcceptanceOptions(
    BasicFlatConverter& ,
    const BasicFlatModelAPI& ma,
    Env& env) {
  auto cal = GetModelAPIAcceptance(ma);
  auto eal = GetModelAPIAcceptanceEXPR(ma);
  auto eial = GetModelAPIAcceptance_EXPR_INTF(ma);
  const bool conacc = (ConstraintAcceptanceLevel::NotAccepted != cal);
  const bool expracc = (ExpressionAcceptanceLevel::NotAccepted != eal);
  const bool expr_intf_acc = (ExpressionAcceptanceLevel::NotAccepted != eial);
  acc_level_item_ = 0;
  if (conacc)
    acc_level_item_
        = std::underlying_type_t<ConstraintAcceptanceLevel>(cal);
  // we prefer expressions, if ModelAPI accepts expression interface
  if (expracc && expr_intf_acc)
    acc_level_item_      // Won't be taken however, if acc:_expr==0
        = std::underlying_type_t<ExpressionAcceptanceLevel>(eal) + 2;
  if (conacc && expracc) {
    env.AddStoredOption(GetAcceptanceOptionNames(),
                        fmt::format(
                            "Solver acceptance level for '{}' as either constraint or expression, "
                            "default {}:\n\n.. value-table::",
                            GetConstraintName(), acc_level_item_).c_str(),
                        acc_level_item_, values_universal_acceptance);
  } else
    if (conacc) {
      env.AddStoredOption(GetAcceptanceOptionNames(),
                          fmt::format(
                              "Solver acceptance level for '{}' as flat constraint, "
                              "default {}:\n\n.. value-table::",
                              GetConstraintName(), acc_level_item_).c_str(),
                          acc_level_item_, values_con_acceptance);
    } else
      if (expracc) {
        env.AddStoredOption(GetAcceptanceOptionNames(),
                            fmt::format(
                                "Solver acceptance level for '{}' as expression, "
                                "default {}:\n\n.. value-table::",
                                GetConstraintName(), acc_level_item_).c_str(),
                            acc_level_item_, values_expr_acceptance);
      } else {
        env.AddStoredOption(GetAcceptanceOptionNames(),
                            "HIDDEN",
                            acc_level_item_, 0, 4);
      }
}

void BasicConstraintKeeper::DoPopulateConstraintList(
    BasicFlatConverter& cvt,
    const BasicFlatModelAPI& ma,
    Env& env) {
  auto cancvt = IfConverterConverts(cvt);
  auto cal = GetModelAPIAcceptance(ma);
  auto eal = GetModelAPIAcceptanceEXPR(ma);
  // Description table
  env.SetConstraintListHeader(
      "List of flat constraints and corresponding expressions.\n"
      "For each constraint/expression, the following are given:\n"
      "\n"
      "  - name,\n"
      "  - convertibility into simpler forms,\n"
      "  - solver acceptance natively as flat constraint,\n"
      "  - solver acceptance natively as expression,\n"
      "  - driver option(s) to modify acceptance\n"
      "    (effective if both convertible and accepted).");
  std::string con_descr = (cancvt) ? "Convertible" : "NonConvertible";
  con_descr += "; ";
  const char * const acc_lev_nam[] = {
       "NotAccepted", "NativeAcceptedButNotRecommended", "NativeRecommended"
  };
  con_descr += acc_lev_nam[std::underlying_type_t<ConstraintAcceptanceLevel>(cal)];
  con_descr += "; ";
  con_descr += acc_lev_nam[std::underlying_type_t<ExpressionAcceptanceLevel>(eal)];
  con_descr += "; ";
  con_descr += GetAcceptanceOptionNames();
  env.AddConstraintDescr(GetConstraintName(), con_descr);
}

template <class Writer>
void WriteVar(Writer& pr, const char* name,
              double lb, double ub, var::Type ty) {
  assert(*name);
  pr << "var " << name;
  if (!lb && 1.0==ub && var::INTEGER==ty)
    pr << " binary";
  else if (lb==ub)
    pr << " = " << lb;
  else {
    if (lb > -DBL_MAX)
    pr << " >=" << lb;
    if (ub < DBL_MAX)
    pr << " <=" << ub;
    if (var::INTEGER == ty)
    pr << " integer";
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const LinTerms& lt,
                    const std::vector<std::string>& vnam) {
  for (int i=0; i<(int)lt.size(); ++i) {
    auto coef = lt.coef(i);
    bool ifpos = coef>=0.0;
    if (i) {
      wrt << (ifpos ? " + " : " - ");
    } else {
      if (!ifpos)
        wrt << "-";
    }
    auto abscoef = std::fabs(coef);
    if (1.0 != abscoef)
      wrt << abscoef << '*';
    wrt << vnam.at(lt.var(i));
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadTerms& qt,
                    const std::vector<std::string>& vnam) {
  for (int i=0; i<(int)qt.size(); ++i) {
    auto coef = qt.coef(i);
    bool ifpos = coef>=0.0;
    if (i) {
      wrt << (qt.coef(i)>=0.0 ? " + " : " - ");
    } else
      if (!ifpos)
        wrt << "-";
    auto abscoef = std::fabs(coef);
    if (1.0 != abscoef)
      wrt << abscoef << '*';
    if (qt.var1(i)==qt.var2(i))
      wrt << vnam.at(qt.var1(i)) << "^2";
    else
      wrt << vnam.at(qt.var1(i))
          << '*' << vnam.at(qt.var2(i));
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadAndLinTerms& qlt,
                    const std::vector<std::string>& vnam) {
  WriteModelItem(wrt, qlt.GetLinTerms(), vnam);
  if (qlt.GetQPTerms().size()) {
    if (qlt.GetLinTerms().size())
      wrt << " + ";
    wrt << '(';
    WriteModelItem(wrt, qlt.GetQPTerms(), vnam);
    wrt << ')';
  }
}

void WriteModelItem(fmt::MemoryWriter& wrt, const QuadraticObjective& obj,
                    const std::vector<std::string>& vnam) {
  wrt << (obj.obj_sense() ? "maximize " : "minimize ");
  assert(obj.name() && *obj.name());
  wrt << obj.name() << ": ";
  WriteModelItem(wrt, obj.GetLinTerms(), vnam);
  if (obj.GetQPTerms().size()) {
    if (obj.GetLinTerms().size())
      wrt << " + ";
    wrt << '(';
    WriteModelItem(wrt, obj.GetQPTerms(), vnam);
    wrt << ')';
  }
}

// Generate
template
void WriteVar(fmt::MemoryWriter& pr, const char* name,
              double lb, double ub, var::Type ty);

template <>
void WriteJSON(JSONW jw, const QuadTerms& qt) {
  jw["coefs"] = qt.coefs();
  jw["vars1"] = qt.vars1();
  jw["vars2"] = qt.vars2();
}

template <>
void WriteJSON(JSONW jw, const LinTerms& qt) {
  jw["coefs"] = qt.coefs();
  jw["vars"] = qt.vars();
}

template <>
void WriteJSON(JSONW jw, const QuadAndLinTerms& qlt) {
  WriteJSON(jw["qp_terms"], qlt.GetQPTerms());
  WriteJSON(jw["lin_terms"], qlt.GetLinTerms());
}

void VisitArguments(const LinTerms& lt, std::function<void (int)> argv) {
  for (auto v: lt.vars())
    argv(v);
}

void VisitArguments(const QuadTerms& lt, std::function<void (int)> argv) {
  for (auto v: lt.vars1())
    argv(v);
  for (auto v: lt.vars2())
    argv(v);
}

void VisitArguments(const QuadAndLinTerms& qlt, std::function<void (int)> argv) {
  VisitArguments(qlt.GetLinTerms(), argv);
  VisitArguments(qlt.GetQPTerms(), argv);
}

/// FlatModelInfo factory
std::unique_ptr<FlatModelInfo> CreateFlatModelInfo() {
  return std::unique_ptr<FlatModelInfo>{new FlatModelInfoImpl()};
}

} // namespace mp
