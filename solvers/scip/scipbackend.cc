#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "scipbackend.h"

extern "C" {
  #include "scip-ampls-c-api.h"    // Scip AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptScip(void* prob) {
  //return SCIP_Interrupt((scip_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateScipBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::ScipBackend()};
}


namespace mp {

/// Create Scip Model Manager
/// @param gc: the Scip common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return ScipModelMgr
std::unique_ptr<BasicModelManager>
CreateScipModelMgr(ScipCommon&, Env&, pre::BasicValuePresolver*&);


ScipBackend::ScipBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateScipModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

ScipBackend::~ScipBackend() {
  CloseSolver();
}

void ScipBackend::OpenSolver() {
  int status = 0;
  // TODO Typically this function creates an instance of the solver environment
  // and an empty model
  void* env_p;
  // Typically try the registered function first;
  // if not available call the solver's API function directly
  /*
  const auto& create_fn = GetCallbacks().cb_initsolver_;
  if (create_fn)
    set_env((GRBenv*)create_fn());
  else
    status = createEnv(&env_p);
    */
  // set_env(env_p);

  /* Todo catch errors
  if ( env() == NULL ) {
    // char  errmsg[CPXMESSAGEBUFSIZE];
    // CPXgeterrorstring (env(), status, errmsg);
     throw std::runtime_error(
       fmt::format("Could not open SCIP environment.\n{}", status) );
  }
  */

  /* TODO Create problem instance
  scip_prob* prob;
  status = SCIP_CreateProb(env_p, &prob);
 */
  Solver::SolverModel* prob = Solver::CreateSolverModel();
  set_lp(prob); // Assign it
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  /* TODO Typically check call */
  /// Turn off verbosity by default
  // SCIP_CCALL(SCIP_SetIntParam(prob, "Logging", 0));

}

void ScipBackend::CloseSolver() {
  /* TODO Cleanup: close problem and environment
  if ( lp() != NULL ) {
    SCIP_CCALL(SCIP_DeleteProb(&lp_) );
  }
  if ( env() != NULL ) {
    SCIP_CCALL(SCIP_DeleteEnv(&env_) );
  }
  */
}

const char* ScipBackend::GetBackendName()
  { return "ScipBackend"; }

std::string ScipBackend::GetSolverVersion() {
  // TODO Return version from solver API
  return "0.0.0";
  //return fmt::format("{}.{}.{}", SCIP_VERSION_MAJOR, 
  //  SCIP_VERSION_MINOR, SCIP_VERSION_TECHNICAL);
}


bool ScipBackend::IsMIP() const {
  // TODO
  return getIntAttr(Solver::VARS_INT) > 0;
  //return getIntAttr(SCIP_INTATTR_ISMIP);
}

bool ScipBackend::IsQCP() const {
  return getIntAttr(Solver::CONS_QUAD) > 0;
// return getIntAttr(SCIP_INTATTR_QELEMS) > 0;
}

ArrayRef<double> ScipBackend::PrimalSolution() {
  int num_vars = NumVars();
  int error;
  std::vector<double> x(num_vars);
  /*
  if (IsMIP()) 
    error = SCIP_GetSolution(lp(), x.data());
  else
    error = SCIP_GetLpSolution(lp(), x.data(), NULL, NULL, NULL);
  if (error)
    x.clear();
    */
  return x;
}

pre::ValueMapDbl ScipBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> ScipBackend::DualSolution_LP() {
  int num_cons = NumLinCons();
  std::vector<double> pi(num_cons);
 // int error = SCIP_GetLpSolution(lp(), NULL, NULL, pi.data(), NULL);
  int error = 0;
  if (error)
    pi.clear();
  return pi;
}

double ScipBackend::ObjectiveValue() const {
 /* if (IsMIP())
    return getDblAttr(SCIP_DBLATTR_BESTOBJ);
  else
    return getDblAttr(SCIP_DBLATTR_LPOBJVAL);
    */
  return 0;
}

double ScipBackend::NodeCount() const {
  return 0;
//  return getIntAttr(SCIP_INTATTR_NODECNT);
}

double ScipBackend::SimplexIterations() const {
  return 0;
//  return getIntAttr(SCIP_INTATTR_SIMPLEXITER);
}

int ScipBackend::BarrierIterations() const {
  return 0;
//  return getIntAttr(SCIP_INTATTR_BARRIERITER);
}

void ScipBackend::ExportModel(const std::string &file) {
  // TODO export proper by file extension
  //SCIP_CCALL(SCIP_WriteLp(lp(), file.data()));
}


void ScipBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptScip, lp());
  // TODO Check interrupter
  //SCIP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void ScipBackend::Solve() {
  if (!storedOptions_.exportFile_.empty()) {
    ExportModel(storedOptions_.exportFile_);
  }
  //SCIP_CCALL(SCIP_Solve(lp()));
  WindupSCIPSolve();
}

void ScipBackend::WindupSCIPSolve() { }

void ScipBackend::ReportResults() {
  ReportSCIPResults();
  BaseBackend::ReportResults();
}

void ScipBackend::ReportSCIPResults() {
  SetStatus( ConvertSCIPStatus() );
  AddSCIPMessages();
  if (need_multiple_solutions())
    ReportSCIPPool();
}
std::vector<double> ScipBackend::getPoolSolution(int i)
{
  std::vector<double> vars(NumVars());
 // SCIP_CCALL(SCIP_GetPoolSolution(lp(), i, NumVars(), NULL, vars.data()));
  return vars;
}
double ScipBackend::getPoolObjective(int i)
{
  double obj;
 // SCIP_CCALL(SCIP_GetPoolObjVal(lp(), i, &obj));
  return obj;
}
void ScipBackend::ReportSCIPPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions;
  /*
  while (++iPoolSolution < getIntAttr(SCIP_INTATTR_POOLSOLS)) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
  */
}


void ScipBackend::AddSCIPMessages() {
  AddToSolverMessage(
          fmt::format("{} simplex iterations\n", SimplexIterations()));
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} barrier iterations\n", nbi));
  if (auto nnd = NodeCount())
    AddToSolverMessage(
          fmt::format("{} branching nodes\n", nnd));
}

std::pair<int, std::string> ScipBackend::ConvertSCIPStatus() {
  namespace sol = mp::sol;
  if (IsMIP())
  {
    /*
    int optstatus = getIntAttr(SCIP_INTATTR_MIPSTATUS);
    switch (optstatus) {
    case SCIP_MIPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case SCIP_MIPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case SCIP_MIPSTATUS_INF_OR_UNB:
      return { sol::INF_OR_UNB, "infeasible or unbounded problem" };
    case SCIP_MIPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case SCIP_MIPSTATUS_TIMEOUT:
    case SCIP_MIPSTATUS_NODELIMIT:
    case SCIP_MIPSTATUS_INTERRUPTED:
      return { sol::INTERRUPTED, "interrupted" };
    }
    */
  }
  else {
    /*
    int optstatus = getIntAttr(SCIP_INTATTR_LPSTATUS);
    switch (optstatus) {
    case SCIP_LPSTATUS_OPTIMAL:
      return { sol::SOLVED, "optimal solution" };
    case SCIP_LPSTATUS_INFEASIBLE:
      return { sol::INFEASIBLE, "infeasible problem" };
    case SCIP_LPSTATUS_UNBOUNDED:
      return { sol::UNBOUNDED, "unbounded problem" };
    case SCIP_LPSTATUS_TIMEOUT:
      return { sol::INTERRUPTED, "interrupted" };
    default:
      return { sol::UNKNOWN, "unfinished" };
    }
    */
  }
  return { sol::UNKNOWN, "not solved" };
}


void ScipBackend::FinishOptionParsing() {
  int v=-1;
 // GetSolverOption(SCIP_INTPARAM_LOGGING, v);
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo lp_values_method[] = {
  { "-1", "Automatic (default)", -1},
  { "1", "Dual simplex", 1},
  { "2", "Barrier", 2},
  { "3", "Crossover", 3},
  { "4", "Concurrent (simplex and barrier simultaneously)", 4},
};


static const mp::OptionValueInfo alg_values_level[] = {
  { "-1", "Automatic (default)", -1},
  { "0", "Off", 0},
  { "1", "Fast", 1},
  { "2", "Normal", 2},
  { "3", "Aggressive", 3}
}; 

static const mp::OptionValueInfo lp_dualprices_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Use Devex pricing algorithm", 0},
  { "1", "Using dual steepest-edge pricing algorithm", 1}
};

static const mp::OptionValueInfo lp_barorder_values_[] = {
  { "-1", "Choose automatically (default)", -1},
  { "0", "Approximate Minimum Degree (AMD)", 0},
  { "1", "Nested Dissection (ND)", 1}
};

void ScipBackend::InitCustomOptions() {

  set_option_header(
      "SCIP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``scip_options``. For example::\n"
      "\n"
      "  ampl: option scip_options 'mipgap=1e-6';\n");

  AddStoredOption("tech:exportfile writeprob writemodel",
      "Specifies the name of a file where to export the model before "
      "solving it. This file name can have extension ``.lp()``, ``.mps``, etc. "
      "Default = \"\" (don't export the model).",
      storedOptions_.exportFile_);

}


double ScipBackend::MIPGap() {
  return 0;
//  return getDblAttr(SCIP_DBLATTR_BESTGAP);
}
double ScipBackend::BestDualBound() {
  return 0;
  //return getDblAttr(SCIP_DBLATTR_BESTBND);
}

double ScipBackend::MIPGapAbs() {
  return std::fabs(
    ObjectiveValue() - BestDualBound());
}


ArrayRef<int> ScipBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  SCIP_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case SCIP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case SCIP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case SCIP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case SCIP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case SCIP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Scip VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> ScipBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  SCIP_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case SCIP_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case SCIP_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case SCIP_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case SCIP_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case SCIP_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown Scip VBasis value: {}", s));
    }
  }*/
  return cons;
}

void ScipBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = SCIP_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = SCIP_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = SCIP_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = SCIP_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = SCIP_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!SCIP_GetColInfo(lp(), SCIP_DBLINFO_LB, 1, index, &lb) && 
        !SCIP_GetColInfo(lp(), SCIP_DBLINFO_UB, 1, index, &ub))
      { 
        if (lb >= -1e-6)
          s = -1;
        else if (ub <= 1e-6)
          s = -2;
        else
          s = -3;  // or, leave at 0?
      }
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL var status value: {}", s));
    }
  }
  SCIP_SetBasis(lp(), stt.data(), NULL);
  */
}

void ScipBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = SCIP_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = SCIP_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  SCIP_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis ScipBackend::GetBasis() {
  std::vector<int> varstt = VarStatii();
  std::vector<int> constt = ConStatii();
  if (varstt.size() && constt.size()) {
    auto mv = GetValuePresolver().PostsolveBasis(
      { std::move(varstt),
        {{{ CG_Linear, std::move(constt) }}} });
    varstt = mv.GetVarValues()();
    constt = mv.GetConValues()();
    assert(varstt.size());
  }
  return { std::move(varstt), std::move(constt) };
}

void ScipBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}


void ScipBackend::ComputeIIS() {
  //SCIP_CCALL(SCIP_ComputeIIS(lp()));
  SetStatus(ConvertSCIPStatus());   // could be new information
}

IIS ScipBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> ScipBackend::VarsIIS() {
  return ArrayRef<int>();
//  return getIIS(lp(), NumVars(), SCIP_GetColLowerIIS, SCIP_GetColUpperIIS);
}
pre::ValueMapInt ScipBackend::ConsIIS() {
  /*auto iis_lincon = getIIS(lp(), NumLinCons(), SCIP_GetRowLowerIIS, SCIP_GetRowUpperIIS);

  std::vector<int> iis_soscon(NumSOSCons());
  SCIP_GetSOSIIS(lp(), NumSOSCons(), NULL, iis_soscon.data());
  ConvertIIS2AMPL(iis_soscon);

  std::vector<int> iis_indicon(NumIndicatorCons());
  SCIP_GetIndicatorIIS(lp(), NumIndicatorCons(), NULL, iis_indicon.data());
  ConvertIIS2AMPL(iis_indicon);

  return { {{ CG_Linear, iis_lincon },
      { CG_SOS, iis_soscon },
      { CG_Logical, iis_indicon }} };
      */
  return { {{ 0, std::vector<int>()}} };
}

void ScipBackend::AddMIPStart(ArrayRef<double> x0) {
  //SCIP_CCALL(SCIP_AddMipStart(lp(), NumVars(), NULL, const_cast<double*>(x0.data())));

}


} // namespace mp


// AMPLs
void* AMPLSOpenScip(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::ScipBackend()},
    slv_opt, cb);
}

void AMPLSCloseScip(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetScipmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::ScipBackend*>(AMPLSGetBackend(slv))->lp();
}
