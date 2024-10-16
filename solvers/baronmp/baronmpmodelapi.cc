#include "baronmpmodelapi.h"

namespace mp {

  std::function<std::string(int)> VExpr::varName = nullptr;

 

  void BaronmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {

    // Have to be written after eventual additional options are injected
    writeBaronOptions();
    VExpr::varName = std::bind(&BaronmpModelAPI::varName, this, std::placeholders::_1);
  }

  void BaronmpModelAPI::AddBounds(const double* bounds, 
    const std::vector<VTYPE>& vtypes,bool lower) {
    bool headerWritten = false;

    double bnd = lower ? -std::numeric_limits<double>::infinity() :
      std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < vtypes.size(); i++)
    {
      if (vtypes[i] == VTYPE::BIN)
        continue;
      if (lower && vtypes[i] == VTYPE::POS)
        continue;
      if (bounds[i] != bnd)
      {
        if (!headerWritten) {
          writeBaron(lower ? "LOWER BOUNDS\{\n" : "UPPER_BOUNDS\{\n");
          headerWritten = true;
        }
        writeBaron(fmt::format("{}: {};\n", lp()->varNames[i], bounds[i]));
      }
    }
    if (headerWritten)
      writeBaron("}\n\n");
  }
  void BaronmpModelAPI::AddVariables(const std::vector<int>& indices, 
    fmt::StringRef prefix, fmt::StringRef header, const char* names) {
    
    if (indices.size() > 0) {
      std::string name;
      writeBaron(fmt::format("{} ", header.to_string()));
      for (int i = 0; i < indices.size(); i++)
      {
        if (i != 0)
          writeBaron(", ");
        if (names)
          name = names[indices[i]];
        else
          name = fmt::format("{}{}", prefix, indices[i]);
        lp()->varMap[name] = indices[i];
        lp()->varNames[indices[i]] = name;
        lp()->baronToAMPLIndices.push_back(indices[i]);

        writeBaron(name);
      }
      writeBaron(";\n");
    }
  }

  void BaronmpModelAPI::AddVariables(const VarArrayDef& v) {
    lp()->varNames.resize(v.size());
    std::vector<int> indicesFree, indicesPos, indicesBin, indicesInt;
    std::vector<VTYPE> vtype(v.size());
    for (auto i = 0; i < v.size(); i++)
    {
      if (v.ptype()[i] == var::Type::INTEGER)
      {
        if ((v.pub()[i] == 1) && (v.plb()[i] == 0))
        {
          vtype[i] = VTYPE::BIN; // BIN
          indicesBin.push_back(i);
        }
        else
        {
          vtype[i] = VTYPE::INT; // INT
          indicesInt.push_back(i);
        }
      }
      else
      {
        if (v.plb()[i] == 0)
        {
          vtype[i] = VTYPE::POS; // POS
          indicesPos.push_back(i);
        }
        else
        {
          vtype[i] = VTYPE::FREE; // FREE
          indicesFree.push_back(i);
        }
      }
    
    }
    int cv = 0;
    std::string name;
    AddVariables(indicesBin, "b", "BINARY_VARIABLES");
    AddVariables(indicesInt, "i", "INTEGER_VARIABLES");
    AddVariables(indicesPos, "p", "POSITIVE_VARIABLES");
    AddVariables(indicesFree, "x", "VARIABLES");
    AddBounds(v.plb(), vtype, true);
    AddBounds(v.pub(), vtype, false);
    nVarsBinary = indicesBin.size();
    nVarsContinuous = indicesPos.size() + indicesFree.size();
    nVarsInteger = indicesInt.size();
  }


  std::string printLinearTerm(double coeff, const char* vname, bool first = false)
  {
    fmt::MemoryWriter w;
    if ((!first) || (coeff > 0))
      w << "+";

    if (std::fabs(coeff) != 1.0) w << coeff << "*";
    w << vname;
    return w.str();
  }
  std::string printQuadTerm(double coeff, const char* v1, const char* v2, bool first = false)
  {
    fmt::MemoryWriter w;
    if ((!first) || (coeff < 0))
      w << " -";

    if (std::fabs(coeff) != 1.0) w << coeff << "*";
    w << v1 << "*" << v2;
    return w.str();
  }


  void BaronmpModelAPI::SetLinearObjective(int iobj, const LinearObjective& lo) {
    if (iobj < 1) {
      fmt::MemoryWriter w;
      w << "OBJ: ";
      if (lo.obj_sense() == mp::obj::Type::MAX)
        w << "maximize ";
      else
        w << "minimize ";

      auto nvars = lo.GetLinTerms().size();
      auto coeffs = lo.GetLinTerms().pcoefs();
      auto vars = lo.GetLinTerms().pvars();

      for (size_t i = 0; i < nvars; i++)
        w << printLinearTerm(coeffs[i], varName(vars[i]).c_str(), i == 0);
      w<<(";\n");
      obj = w.str();
     
    }
    else {
      // TODO
      fmt::print("Setting {}-th linear objective: {} terms.\n", iobj, lo.num_terms());
    }
  }


  std::string BaronmpModelAPI::varName(int index) {
    return lp()->varNames[index];
  }

  std::string BaronmpModelAPI::createConName(const std::string& name) {
    std::string n;
    if (name.empty())
      n = fmt::format("c{}", ++n_unamed_constraint);
    else
      n = name;
    lp()->conNames.push_back(n);
    return n;
  }
  
  std::string addRhs(const std::string& body, double lhs, double rhs) {
    fmt::MemoryWriter w;
    bool hasLhs = lhs != -std::numeric_limits<double>::infinity();
    bool hasRhs = rhs != std::numeric_limits<double>::infinity();
    if ((hasLhs && hasRhs) && (lhs != rhs))
      w << fmt::format("{} <= ", lhs);
    w << body;
    if (hasRhs)
      w << fmt::format(" {} {}", (lhs != rhs) ? "<=" : "==", rhs);
    else if (hasLhs)
      w << fmt::format(" >= {}", lhs);
    w << ";\n";
    return w.str();

  }
  void BaronmpModelAPI::addLinear(const std::string& name,
    double lhs, double rhs,
    size_t nvars, const int* vars, const double* coeffs){
    

    fmt::MemoryWriter body;
    for (size_t i = 0; i < nvars; i++)
      body << printLinearTerm(coeffs[i], varName(vars[i]).c_str(), i == 0);
    fmt::MemoryWriter w;
    w << createConName(name) << ": ";
    w << addRhs(body.str(), lhs, rhs);
    cons.push_back(w.str());
}

  void BaronmpModelAPI::addQuadratic(const std::string& name,
    double lhs, double rhs, const mp::LinTerms& lt,
    const mp::QuadTerms& qt) {
    fmt::MemoryWriter body;
    for (size_t i = 0; i < lt.size(); i++)
      body << printLinearTerm(lt.pcoefs()[i], varName(lt.var(i)).c_str(), i == 0);
    for (size_t i = 0; i < qt.size(); i++)
      body << printQuadTerm(qt.coef(i), varName(qt.var1(i)).c_str(), varName(qt.var2(i)).c_str(), i == 0);
    fmt::MemoryWriter w;
    w << createConName(name) << ": ";
    w << addRhs(body.str(), lhs, rhs);
    cons.push_back(w.str());
  }

  template <typename Args, typename Params, typename NumOrLogic, typename Id>
  void BaronmpModelAPI::addFunctionalConstraint(const std::string& func,
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c) {
    fmt::MemoryWriter w;
    
    w<< createConName(c.GetName()) << ": ";

    w << fmt::format("{} = {}({}", varName(c.GetResultVar()), func, varName(c.GetArguments()[0]));
    if (c.GetArguments().size() > 1) {
      for (auto i = 1; i < c.GetArguments().size() - 1; ++i)
        w << fmt::format("{}, ", varName(c.GetArguments()[i]));
      w << varName(c.GetArguments()[c.GetArguments().size() - 1]);
    }
    w << ");\n";
    cons.push_back(w.str());

  }


void BaronmpModelAPI::AddConstraint(const LinConRange& lc) {
  addLinear(lc.name(), lc.lb(), lc.ub(),
    lc.size(), lc.pvars(), lc.pcoefs()) ;
}
void BaronmpModelAPI::AddConstraint(const LinConLE& lc) {
  addLinear(lc.name(), 
    -std::numeric_limits<double>::max(),lc.ub(),
    lc.size(), lc.pvars(), lc.pcoefs());
}
void BaronmpModelAPI::AddConstraint(const LinConEQ& lc) {
  addLinear(lc.name(),
    lc.lb(), lc.ub(),
    lc.size(), lc.pvars(), lc.pcoefs());
}
void BaronmpModelAPI::AddConstraint(const LinConGE& lc) {
  addLinear(lc.name(),
    lc.lb(), std::numeric_limits<double>().max(),
    lc.size(), lc.pvars(), lc.pcoefs());
}

void BaronmpModelAPI::AddConstraint(const QuadConRange& qc) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint( const QuadConLE& qc ) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint( const QuadConEQ& qc ) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint( const QuadConGE& qc ) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint(const MaxConstraint& mc) {
  addFunctionalConstraint("max", mc);
}

void BaronmpModelAPI::AddConstraint(const MinConstraint& mc) {
  addFunctionalConstraint("min", mc);
}

void BaronmpModelAPI::AddConstraint(const AbsConstraint& c) {
  assert(c.GetArguments().size() == 1);
  fmt::MemoryWriter w;
  w << createConName(c.GetName()) << ": ";
  w << fmt::format("{} = (({}", varName(c.GetResultVar()), varName(c.GetArguments()[0]));
  w << "^2 ) ^ 0.5);\n";
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const AndConstraint& cc) {
  addFunctionalConstraint("and", cc);
}
void BaronmpModelAPI::AddConstraint(const OrConstraint& dc) {
  addFunctionalConstraint("or", dc);
}

void BaronmpModelAPI::AddConstraint(const ExpConstraint& cc) {
  addFunctionalConstraint("exp", cc);
}

void BaronmpModelAPI::AddConstraint(const LogConstraint& cc) {
  addFunctionalConstraint("log", cc);
}

void BaronmpModelAPI::AddConstraint(const ExpAConstraint& cc) {
  // TODO
  fmt::print("Adding ExpA constraint \"{}\"\n", cc.GetName());
}



void BaronmpModelAPI::AddConstraint(const LogAConstraint& cc) {
  // TODO
  fmt::print("Adding LogA constraint \"{}\"\n", cc.GetName());
}

void BaronmpModelAPI::AddConstraint(const PowConstraint& cc) {
  // TODO
}



void BaronmpModelAPI::AddConstraint(const PLConstraint& plc) {
  fmt::print("Adding PL constraint \"{}\"\n", plc.GetName());
  // todo
}


void BaronmpModelAPI::FinishProblemModificationPhase() {
  if (lp()->conNames.size() != 0)
  {
    writeBaron("EQUATIONS ");
    writeBaron(lp()->conNames[0]);
    for (auto i = 0; i < lp()->conNames.size(); ++i)
      writeBaron(" ,{}", lp()->conNames[i]);
    writeBaron(";\n\n");

    for (auto c : cons)
      writeBaron(c);
  }
  writeBaron("\n");
  writeBaron(obj);
}


template <int SENSE> void BaronmpModelAPI::addTopLevel(const NLBaseAssign<SENSE>& c) {
  const auto var = GetVariable(c);
  fmt::MemoryWriter w;
  w << createConName(c.GetName()) << ": ";
  const std::string sense[] = { "<=", "=", ">="};
  w << fmt::format("{} {} ", varName(var), sense[SENSE+1]);
  const auto& frm = GetExpression(c);
  w << frm.ToString(false); 
  w << ";\n";
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const NLConstraint& nlc) {
  const auto& exp = GetExpression(nlc);
  fmt::MemoryWriter body;
  for (int i = 0; i < GetLinSize(nlc); ++i)
    body << printLinearTerm(GetLinCoef(nlc, i), varName(GetLinVar(nlc, i)).c_str(), i == 0);
  if (GetLinSize(nlc) > 0) body <<" + ";
  body << exp.ToString(false);
  fmt::MemoryWriter w;
  w << createConName(nlc.GetName()) << ": ";
  w << addRhs(body.str(), GetLower(nlc), GetUpper(nlc));
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const NLAssignEQ& nla) {
  addTopLevel(nla);
}
void BaronmpModelAPI::AddConstraint(const NLAssignLE& nla) {
  addTopLevel(nla);
}
void BaronmpModelAPI::AddConstraint(const NLAssignGE& nla) {    
  addTopLevel(nla);
}

VExpr BaronmpModelAPI::AddExpression(const NLAffineExpression& nla) {
  VExpr tree;
  tree.opcode = Opcode::ADD;
  AppendLinAndConstTerms(tree, nla);
  return tree;
}

VExpr  BaronmpModelAPI::AddExpression(const NLQuadExpression& nlq) {
  VExpr tree;
  tree.opcode = Opcode::ADD;
  AppendLinAndConstTerms(tree, nlq);
  AppendQuadTerms(tree, nlq);
  return tree;
}

template <class MPExpr>
void BaronmpModelAPI::AppendLinAndConstTerms(
  Expr &ff, const MPExpr& nla) {
  if (double ct = GetConstTerm(nla))
    ff.AddArgument(VExpr::makeConstant(ct));

  for (int i = 0; i < GetLinSize(nla); ++i) {
    if (double coef = GetLinCoef(nla, i)) {
      if (1.0 != coef) {
        auto f1 = VExpr(Opcode::MULTIPLY);
        f1.AddArgument(VExpr::makeConstant(coef));
        f1.AddArgument(GetLinTerm(nla, i));
        ff.AddArgument(f1);
      }
      else {
        ff.AddArgument(GetLinTerm(nla, i));
      }
    }
  }
}


template <class MPExpr>
void BaronmpModelAPI::AppendQuadTerms(
  VExpr& ff, const MPExpr& nlq) {
  for (int i = 0; i < GetQuadSize(nlq); ++i) {
    if (double coef = GetQuadCoef(nlq, i)) {
      auto f1 = VExpr(Opcode::MULTIPLY);
      f1.AddArgument(GetQuadTerm1(nlq, i));
      f1.AddArgument(GetQuadTerm2(nlq, i));
      if (1.0 != coef) {
        f1.AddArgument(VExpr::makeConstant(coef));
      }
      ff.AddArgument(f1);
    }
  }
}
VExpr BaronmpModelAPI::AddExpression(const ExpAExpression& e) {
  auto res = VExpr(Opcode::POW);
  // Constant^Variable
  res.AddArgument(VExpr::makeConstant(GetParameter(e, 0)));
  res.AddArgument(GetArgExpression(e, 0));
  return res;
}

VExpr BaronmpModelAPI::AddExpression(const ExpExpression& e) {
  // e^v
  return VExpr(Opcode::EXP, GetArgExpression(e, 0));
}

VExpr BaronmpModelAPI::AddExpression(const LogExpression& e) {
  return VExpr(Opcode::LOG, GetArgExpression(e, 0));
}
VExpr BaronmpModelAPI::AddExpression(const LogAExpression& e) {
  auto logvar = VExpr(Opcode::LOG, GetArgExpression(e, 0));
  auto logbase = VExpr(Opcode::LOG, VExpr::makeConstant(GetParameter(e, 0)));
  auto expr = VExpr(Opcode::DIVIDE);
  expr.AddArgument(logvar);
  expr.AddArgument(logbase);
  return expr;
}

VExpr BaronmpModelAPI::AddExpression(const PowExpression& e) {
  auto expr=VExpr(Opcode::POW, GetArgExpression(e, 0));
  for (int i = 0; i < GetNumParameters(e); ++i)
    expr.AddArgument(VExpr::makeConstant(GetParameter(e, i)));
  return expr;
}

VExpr BaronmpModelAPI::AddExpression(const AbsExpression& e) {
  auto f1 = VExpr(Opcode::POW, GetArgExpression(e, 0)); // x^2
  f1.AddArgument(VExpr::makeConstant(2));
  auto expr = VExpr(Opcode::POW, f1); // (x^2)^0.5
  expr.AddArgument(VExpr::makeConstant(0.5));
  return expr;
}






} // namespace mp
