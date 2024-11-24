#include "cbcmpmodelapi.h"

#include <algorithm> // for std::copy

namespace mp {

void CbcmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {
}

std::string sanitizeName(std::string n) {
  // Cbc does not like square brackets or spaces in variable names
  std::replace(n.begin(), n.end(), '[', '(');
  std::replace(n.begin(), n.end(), ']', ')');
  std::replace(n.begin(), n.end(), ' ', '_');
  return n;
}

void CbcmpModelAPI::AddVariables(const VarArrayDef& v) {
  m.SetNumVars(v.size());
  std::copy(v.plb(), v.plb() + v.size(), m.lb.begin());
  std::copy(v.pub(), v.pub() + v.size(), m.ub.begin());
  for (int i = 0; i < v.size(); i++) {
    m.varnames[i] = v.pnames() ? sanitizeName(v.pnames()[i]) : fmt::format("x{}", i);
    if (var::Type::INTEGER == v.ptype()[i])
      m.varinteger.push_back(i);
  }
}

void CbcmpModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  m.objsense = lo.obj_sense() == obj::MAX ? -1 : 1;
  if (iobj<1) {
    for (int i = 0; i < lo.num_terms(); i++)
      m.obj[lo.vars()[i]] = lo.coefs()[i];
  } else {
    throw std::runtime_error("Multiple objectives not supported");
  }
}

void CbcmpModelAPI::AddConstraint(const LinConLE& lc) {
  m.AddRow(lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), -Infinity(), lc.rhs());
}
void CbcmpModelAPI::AddConstraint(const LinConEQ& lc) {
  m.AddRow(lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), lc.rhs(), lc.rhs());
}
void CbcmpModelAPI::AddConstraint(const LinConGE& lc) {
  m.AddRow(lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), lc.rhs(), Infinity());
}
void CbcmpModelAPI::AddConstraint(const LinConRange& lc) {
  m.AddRow(lc.GetName(), lc.size(), lc.pvars(),
    lc.pcoefs(), lc.lb(), lc.ub());
}

void CbcmpModelAPI::AddConstraint(const SOS1Constraint& sos) {
  m.addSOS(sos.get_vars().size(),
    sos.get_vars().data(), sos.get_weights().data(), 1);
}

void CbcmpModelAPI::AddConstraint(const SOS2Constraint& sos) {
  m.addSOS(sos.get_vars().size(),
    sos.get_vars().data(), sos.get_weights().data(), 2);

}

void CbcmpModelAPI::FinishProblemModificationPhase() {
  // Push all model to CBC
  std::vector<double> values;
  std::vector<int> row_indices;
  std::vector<CoinBigIndex> col_pointers;
  m.matrix.toCSC(values, row_indices, col_pointers);
  Cbc_loadProblem(lp(),
    m.lb.size(), m.lhs.size(),
    col_pointers.data(),
    row_indices.data(), values.data(), m.lb.data(), m.ub.data(),
    m.obj.data(), m.lhs.data(), m.rhs.data());
  Cbc_setObjSense(lp(), m.objsense);
  for(auto i : m.varinteger) Cbc_setInteger(lp(), i);
  // If SOS constraints are present, push those also
  for(int i=0; i<2; i++) {
    if (m.soslist[i].sos_rowstarts.size() > 0) {
      m.soslist[i].sos_rowstarts.push_back(m.soslist[i].sos_colindices.size());
      Cbc_addSOS(lp(), m.soslist[i].sos_types.size(),
        m.soslist[i].sos_rowstarts.data(), m.soslist[i].sos_colindices.data(),
        m.soslist[i].sos_colweights.data(), i+1);

    }
  }
  for (auto i = 0; i < m.varnames.size(); i++) 
    if(!m.varnames[i].empty())
        Cbc_setColName(lp(), i, m.varnames[i].c_str());

  for (auto i = 0; i < m.connames.size(); i++)
    if (!m.connames[i].empty())
      Cbc_setRowName(lp(), i, m.connames[i].c_str());
  m.Clear(); // release memory before solving

}


} // namespace mp
