#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "mp2nlbackend.h"

extern "C" {
  #include "mp2nl-ampls-c-api.h"    // MP2NL AMPLS C API
}
#include "mp/ampls-cpp-api.h"


std::unique_ptr<mp::BasicBackend> CreateMP2NLBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::MP2NLBackend()};
}


namespace mp {

/// Create MP2NL Model Manager
/// @param gc: the MP2NL common handle
/// @param e: environment
/// @param pre: presolver to be returned,
/// need it to convert solution data
/// @return MP2NLModelMgr
std::unique_ptr<BasicModelManager>
CreateMP2NLModelMgr(MP2NLCommon&, Env&, pre::BasicValuePresolver*&);


MP2NLBackend::MP2NLBackend() {
  OpenSolver();

  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateMP2NLModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);

  p_nlsi_ = get_other()->p_nlsi_;    // before copying
  /// Copy env/lp to ModelAPI
  copy_common_info_to_other();
}

MP2NLBackend::~MP2NLBackend() {
  CloseSolver();
}


const char* MP2NLBackend::GetBackendName()
  { return "MP2NLBackend"; }

std::string MP2NLBackend::GetSolverVersion() {
  return "0.1";
}

std::string MP2NLBackend::set_external_libs() {
  return "NLWriter2";
}

ArrayRef<double> MP2NLBackend::PrimalSolution() {
  return GetNLSolver().GetX();
}

pre::ValueMapDbl MP2NLBackend::DualSolution() {
  return {{ { CG_Algebraic, DualSolution_LP() } }};
}

ArrayRef<double> MP2NLBackend::DualSolution_LP() {
  return GetNLSolver().GetY();
}

double MP2NLBackend::ObjectiveValue() const {
  return 0.0;                        // SOL does not provide one
}

double MP2NLBackend::NodeCount() const {
  return 0.0;
}

double MP2NLBackend::SimplexIterations() const {
  return 0.0;
}

int MP2NLBackend::BarrierIterations() const {
  return 0;
}

void MP2NLBackend::ExportModel(const std::string &file) {
}


void MP2NLBackend::SetInterrupter(mp::Interrupter *inter) {
  // TODO
}

void MP2NLBackend::Solve() {
  GetNLSolver().Solve(
      storedOptions_.solver_.c_str(),
      storedOptions_.solver_options_.c_str());
  WindupMP2NLSolve();
}

void MP2NLBackend::WindupMP2NLSolve() { }

void MP2NLBackend::ReportResults() {
  ReportMP2NLResults();
  BaseBackend::ReportResults();
}

void MP2NLBackend::ReportMP2NLResults() {
  SetStatus( GetSolveResult() );
  AddMP2NLMessages();
  if (need_multiple_solutions())
    ReportMP2NLPool();
}
std::vector<double> MP2NLBackend::getPoolSolution(int i)
{
  int num_vars = NumVars();
  std::vector<double> vars(num_vars);
  return vars;
}
double MP2NLBackend::getPoolObjective(int i)
{
  double obj {0.0};
  return obj;
}
void MP2NLBackend::ReportMP2NLPool() {
  if (!IsMIP())
    return;
  int iPoolSolution = -1;
  int nsolutions = 0;
  
  while (++iPoolSolution < nsolutions) {
    ReportIntermediateSolution(
      { getPoolSolution(iPoolSolution),
        {}, { getPoolObjective(iPoolSolution) } });
  }
}


void MP2NLBackend::AddMP2NLMessages() {
}

std::pair<int, std::string> MP2NLBackend::GetSolveResult() {
  return {
          GetNLSolver().GetSolveResult(),
      GetNLSolver().GetSolveMessage() };
}


void MP2NLBackend::FinishOptionParsing() {
  int v=storedOptions_.outlev_;
  set_verbose_mode(v>0);
}


////////////////////////////// OPTIONS /////////////////////////////////

static const mp::OptionValueInfo childsel[] = {
  { "d", "down", 0},
  { "u", "up", 1},
  { "p", "pseudo costs", 2},
  { "i", "inference", 3},
  { "l", "lp value", 4},
  { "r", "root LP value difference", 5},
  { "h", "hybrid inference/root LP value difference (default)", 6}
};


void MP2NLBackend::InitCustomOptions() {

  set_option_header(
    "MP2NL Optimizer Options for AMPL\n"
    "--------------------------------------------\n"
    "\n"
    "To set these options, assign a string specifying their values to the "
    "AMPL option ``mp2nl_options``. For example::\n"
    "\n"
    "  ampl: option mp2nl_options 'solver=baron solver_options=\"outlev=1 iisfind=1\"';\n");

  AddStoredOption("tech:outlev outlev",
    "0*/1: Verbosity for the MP2NL driver. For the underlying solver, use tech:solver_options.",
    storedOptions_.outlev_);

  // AddStoredOption("tech:logfile logfile",
  //                 "Log file name.",
  //                 storedOptions_.logFile_);

  AddStoredOption("tech:solver solver",
                  "Underlying AMPL solver.",
                  storedOptions_.solver_);

  AddStoredOption("tech:solver_options solver_options slv_opts",
                  "Underlying solver options.",
                  storedOptions_.solver_options_);

  /// Enforce time limit?
  AddStoredOption("lim:time timelim timelimit time_limit",
    "Enforce limit on solve time (in seconds; default: 1e+20). DOES NOTHING CURRENTLY.",
    storedOptions_.tilim_);


  /// Custom infinity value?
  // AddSolverOption("num:infinity infinity",
  //   "Values larger than this are considered infinity (default: 1e+20)",
  //   "numerics/infinity", 1e+10, DBL_MAX);

}


double MP2NLBackend::MIPGap() {
  return 0.0;
}
double MP2NLBackend::BestDualBound() {
  return 0.0;
}

double MP2NLBackend::MIPGapAbs() {
  return 0.0;
}


ArrayRef<int> MP2NLBackend::VarStatii() {
  
  std::vector<int> vars(NumVars());
  /*
  MP2NL_GetBasis(lp(), vars.data(), NULL);
  for (auto& s : vars) {
    switch (s) {
    case MP2NL_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case MP2NL_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case MP2NL_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case MP2NL_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case MP2NL_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown MP2NL VBasis value: {}", s));
    }
  }
  */
  return vars;
}

ArrayRef<int> MP2NLBackend::ConStatii() {

  std::vector<int> cons(NumLinCons());
  /*
  MP2NL_GetBasis(lp(), NULL, cons.data());
  for (auto& s : cons) {
    switch (s) {
    case MP2NL_BASIS_BASIC:
      s = (int)BasicStatus::bas;
      break;
    case MP2NL_BASIS_LOWER:
      s = (int)BasicStatus::low;
      break;
    case MP2NL_BASIS_UPPER:
      s = (int)BasicStatus::upp;
      break;
    case MP2NL_BASIS_SUPERBASIC:
      s = (int)BasicStatus::sup;
      break;
    case MP2NL_BASIS_FIXED:
      s = (int)BasicStatus::equ;
      break;
    default:
      MP_RAISE(fmt::format("Unknown MP2NL VBasis value: {}", s));
    }
  }*/
  return cons;
}

void MP2NLBackend::VarStatii(ArrayRef<int> vst) {
  int index[1];
  std::vector<int> stt(vst.data(), vst.data() + vst.size());
  /*
  for (auto j = stt.size(); j--; ) {
    auto& s = stt[j];
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = MP2NL_BASIS_BASIC;
      break;
    case BasicStatus::low:
      s = MP2NL_BASIS_LOWER;
      break;
    case BasicStatus::equ:
      s = MP2NL_BASIS_FIXED;
      break;
    case BasicStatus::upp:
      s = MP2NL_BASIS_UPPER;
      break;
    case BasicStatus::sup:
    case BasicStatus::btw:
      s = MP2NL_BASIS_SUPERBASIC;
      break;
    case BasicStatus::none:
      /// 'none' is assigned to new variables. Compute low/upp/sup:
      /// Depending on where 0.0 is between bounds
      double lb, ub;
      index[0] = (int)j;
      if(!MP2NL_GetColInfo(lp(), MP2NL_DBLINFO_LB, 1, index, &lb) &&
        !MP2NL_GetColInfo(lp(), MP2NL_DBLINFO_UB, 1, index, &ub))
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
  MP2NL_SetBasis(lp(), stt.data(), NULL);
  */
}

void MP2NLBackend::ConStatii(ArrayRef<int> cst) {
  /*
  std::vector<int> stt(cst.data(), cst.data() + cst.size());
  for (auto& s : stt) {
    switch ((BasicStatus)s) {
    case BasicStatus::bas:
      s = MP2NL_BASIS_BASIC;
      break;
    case BasicStatus::none:   // for 'none', which is the status
    case BasicStatus::upp:    // assigned to new rows, it seems good to guess
    case BasicStatus::sup:    // a valid status.
    case BasicStatus::low:    // 
    case BasicStatus::equ:    // For active constraints, it is usually 'sup'.
    case BasicStatus::btw:    // We could compute slack to decide though.
      s = MP2NL_BASIS_SUPERBASIC;
      break;
    default:
      MP_RAISE(fmt::format("Unknown AMPL con status value: {}", s));
    }
  }
  MP2NL_SetBasis(lp(), NULL, stt.data());
  */
}

SolutionBasis MP2NLBackend::GetBasis() {
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

void MP2NLBackend::SetBasis(SolutionBasis basis) {
  auto mv = GetValuePresolver().PresolveBasis(
    { basis.varstt, basis.constt });
  auto varstt = mv.GetVarValues()();
  auto constt = mv.GetConValues()(CG_Linear);
  assert(varstt.size());
  assert(constt.size());
  VarStatii(varstt);
  ConStatii(constt);
}

void MP2NLBackend::AddPrimalDualStart(Solution sol)
{
  auto mv = GetValuePresolver().PresolveSolution(
        { sol.primal, sol.dual } );
  auto x0 = mv.GetVarValues()();
	auto pi0 = mv.GetConValues()(CG_Linear);
}

void MP2NLBackend::AddMIPStart(ArrayRef<double> x0, ArrayRef<int> sparsity) {
}

void MP2NLBackend::DoWriteProblem(const std::string& name) {
  ExportModel(name);
}

} // namespace mp


// AMPLs
AMPLS_MP_Solver* Open_MP2NL(CCallbacks cb = {}) {
  return AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend>{new mp::MP2NLBackend()},
    cb);
}

void AMPLSClose_MP2NL(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* AMPLSGetModel_MP2NL(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::MP2NLBackend*>(AMPLSGetBackend(slv));
}
