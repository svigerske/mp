#include <cmath>

#include "highsmpmodelapi.h"


namespace mp {

void HighsModelAPI::InitProblemModificationPhase(const FlatModelInfo*) { }

void HighsModelAPI::AddVariables(const VarArrayDef& v) {
  std::vector<int> intIndices;
  std::vector<double> costs(v.size(), 0);
  std::vector<double> lbs(v.size());
  std::vector<double> ubs(v.size());
  for (int i = 0; i < v.size(); i++) {
    if (var::Type::INTEGER == v.ptype()[i]) intIndices.push_back(i);
    lbs[i] = std::isinf(v.plb()[i]) ? MinusInfinity() : v.plb()[i];
    ubs[i] = std::isinf(v.pub()[i]) ?  Infinity() : v.pub()[i];

  }
  HIGHS_CCALL(Highs_addCols(lp(), v.size(), costs.data(), lbs.data(), ubs.data(), 0, NULL, NULL, NULL));
  if (intIndices.size() > 0) {
    std::vector<int> types(intIndices.size(), 1); // TODO get the 1 from solver API?
    HIGHS_CCALL(Highs_changeColsIntegralityBySet(lp(), intIndices.size(),
      intIndices.data(), types.data()));
  }
  if (v.pnames())
    for (int i = 0; i < v.size(); i++)
      HIGHS_CCALL(Highs_passColName(lp(), i, v.pnames()[i]));
}

void HighsModelAPI::SetLinearObjective( int iobj, const LinearObjective& lo ) {
  if (iobj < 1) {
    // Highs_changeObjectiveOffset(highs, offset); TODO offset?
    std::vector<double> objc_new(NumVars());         // dense vector to
    for (auto i=lo.vars().size(); i--; )
      objc_new[lo.vars()[i]] = lo.coefs()[i];
    HIGHS_CCALL(Highs_changeColsCostByRange(lp(), 0, NumVars()-1, objc_new.data()));
    HIGHS_CCALL(Highs_changeObjectiveSense(lp(), 
                                           obj::Type::MAX==lo.obj_sense() ?
                                               kHighsObjSenseMaximize : kHighsObjSenseMinimize) );
  } else {
      throw std::runtime_error("HighS does not support multiple objectives.");
  }
}


void HighsModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
  if (1 > iobj) {
    SetLinearObjective(iobj, qo);
    const auto& qt = qo.GetQPTerms();
    std::vector<int> startCols(NumVars());
    std::vector<double> coeffs(qt.size());
    // Convert to Highs Hessian upper triangular format.
    size_t q=0;              // the index in qt
    for (size_t j=0; j<startCols.size(); ++j) {
      assert(q>=qt.size() || j<=(size_t)qt.var1(q)); // qt sorted
      if (q<qt.size() && j==(size_t)qt.var1(q)) {
        startCols[j] = q;
        for ( ; q<qt.size() && (size_t)qt.var1(q) == j; ++q) {
          assert(j <= (size_t)qt.var2(q));      // upper triangular
          coeffs[q] =
                (j == (size_t)qt.var2(q) ? 2.0 : 1.0) * qt.coef(q);
        }
      } else {
        startCols[j] = q;
        // Could do these but not needed:
        //        index.push_back(j);
        //        coeffs.push_back(0.0);  // empty diagonal element
      }
    }
    HIGHS_CCALL(Highs_passHessian(lp(), NumVars(), qt.size(),
                                  kHighsHessianFormatTriangular,
      startCols.data(), qt.pvars2(), coeffs.data()));
  }
  else {
    throw std::runtime_error("Multiple quadratic objectives not supported");
  }
}

void HighsModelAPI::AddConstraint(const LinConRange& lc) {
  acc_constraints_.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConLE& lc) {
  acc_constraints_.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConEQ& lc) {
  acc_constraints_.add(lc);
}
void HighsModelAPI::AddConstraint(const LinConGE& lc) {
  acc_constraints_.add(lc);
}

void HighsModelAPI::FinishProblemModificationPhase() {
  HIGHS_CCALL(Highs_addRows(lp(),
    acc_constraints_.lb.size(),
    acc_constraints_.lb.data(),
    acc_constraints_.ub.data(),
    acc_constraints_.coeffs.size(),
    acc_constraints_.starts.data(),
    acc_constraints_.indices.data(),
    acc_constraints_.coeffs.data()));
  acc_constraints_ = AccConstraints();      // reinitialize accumulator for model modification
}

} // namespace mp
