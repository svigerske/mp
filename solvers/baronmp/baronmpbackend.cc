#include <vector>
#include <climits>
#include <cfloat>

#include "mp/env.h"
#include "mp/flat/model_api_base.h"
#include "baronmpbackend.h"

extern "C" {
  #include "baronmp-ampls-c-api.h"    // Baronmp AMPLS C API
}
#include "mp/ampls-cpp-api.h"

namespace {


bool InterruptBaronmp(void* prob) {

  //return BARONMP_Interrupt((baronmp_prob*)prob);
  return true;
}

}  // namespace {}

std::unique_ptr<mp::BasicBackend> CreateBaronmpBackend() {
  return std::unique_ptr<mp::BasicBackend>{new mp::BaronmpBackend()};
}


namespace mp {

/// Create Baronmp Model Manager
/// @param gc: the Baronmp common handle
/// @param e: environment
/// @param pre: presolver to be names[Solver::returned]= "";
/// need it to convert solution data
/// @return BaronmpModelMgr
std::unique_ptr<BasicModelManager>
CreateBaronmpModelMgr(BaronmpCommon&, Env&, pre::BasicValuePresolver*&);


BaronmpBackend::BaronmpBackend() {
  OpenSolver();
  /// Create a ModelManager
  pre::BasicValuePresolver* pPre;
  auto data = CreateBaronmpModelMgr(*this, *this, pPre);
  SetMM( std::move( data ) );
  SetValuePresolver(pPre);



}

BaronmpBackend::~BaronmpBackend() {
  CloseSolver();
}

void BaronmpBackend::OpenSolver() {
  set_lp(new BaronProblemInfo());
  int status = 0;
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
}

void BaronmpBackend::CloseSolver() {
  if(!baronOptions().keepsol)
  recrmdir(baronDir);
} 

const char* BaronmpBackend::GetBackendName()
  { return "BaronmpBackend"; }

std::string BaronmpBackend::GetSolverVersion() {
  return fmt::format("{}.{}.{}", v_year, v_month, v_day);
}


bool BaronmpBackend::IsMIP() const {
  return false;
}

bool BaronmpBackend::IsQCP() const {
  return false;
}

ArrayRef<double> BaronmpBackend::PrimalSolution() {
  auto vec = resFileData_.PrimalSolution();
  int num_vars = vec.size();
  int error;
  std::vector<double> x(num_vars);
  for (int i = 0; i < num_vars; i++)
    x[lp()->baronToAMPLIndices[i]] = vec[i];
  return x;
}

pre::ValueMapDbl BaronmpBackend::DualSolution() {
  return {{ { CG_Linear, DualSolution_LP() } }};
}

ArrayRef<double> BaronmpBackend::DualSolution_LP() {
  return resFileData_.DualSolution();
}

double BaronmpBackend::ObjectiveValue() const {
  return resFileData_.ObjectiveValue();
}

double BaronmpBackend::NodeCount() const {
  return 0;
}

double BaronmpBackend::SimplexIterations() const {
  return 0;
}

int BaronmpBackend::BarrierIterations() const {
  return timFileData_.itera;
}


void BaronmpBackend::SetInterrupter(mp::Interrupter *inter) {
  inter->SetHandler(InterruptBaronmp, lp());
  // TODO Check interrupter
  //BARONMP_CCALL( CPXsetterminate (env(), &terminate_flag) );
}

void BaronmpBackend::Solve() {
  FILE_BAR->close();

  GlobalsManager::BaronGlobals g;
  g.baronDir = baronDir;
  g.initialDir = initialDir;
  g.lsolmsg = baronOptions().lsolmsg;
  if (!baronOptions().lsolver.empty())
    g.lsolver = baronOptions().lsolver;
  g.lpsol_pvsub = baronOptions().lpsol;
  if (baronOptions().lpsolver == "cplex")
    g.lpsol_dll = baronOptions().cplexlibname;
  g.nlfile = nlFilePath;
  g.nvars = nVarsInteger + nVarsContinuous + nVarsBinary;
  g.verbuf = version();
  g.v_keepsol = baronOptions().keepsol;
  GlobalsManager::writeGlobals(FILENAME_AMPL, g);
  int status =  BaronmpCommon::runBaron(filePathBar);
  if (status)
    errorLevel = -1; // Problems running process

  WindupBARONMPSolve();
}

void BaronmpBackend::WindupBARONMPSolve() {
  timFileData_ = TimFileData::fromFile(BARON_TIM);
  resFileData_ = ResFileData::fromFile(BARON_RES, timFileData_.nodeopt);

}

void BaronmpBackend::ReportResults() {
  ReportBARONMPResults();
  BaseBackend::ReportResults();
}

void BaronmpBackend::ReportBARONMPResults() {
  SetStatus( GetSolveResult() );
  AddBARONMPMessages();
  if (need_multiple_solutions())
    ReportBARONMPPool();
}

void BaronmpBackend::ReportBARONMPPool() {
  for(int i=0; i<resFileData_.NumSol();  i++)
  for(const auto &sol : resFileData_.PrimalSolutions()){
    ReportIntermediateSolution(
      { resFileData_.PrimalSolutions()[i],
        {}, {  resFileData_.Objectives()[i]}});
      }
}

void BaronmpBackend::AddBARONMPMessages() {
  if (auto nbi = BarrierIterations())
    AddToSolverMessage(
          fmt::format("{} iterations\n", nbi));
}

std::pair<int, std::string> BaronmpBackend::GetSolveResult() {
  namespace sol = mp::sol;

  if(errorLevel == -1) // could not call Baron process
    return { sol::UNKNOWN, "not solved, could not call Baron process" };
  switch (timFileData_.modelstatus) {
  case 1:
    return { sol::SOLVED, "optimal solution" };
  case 2:
    return { sol::INFEASIBLE, "infeasible problem" };
  case 3:
    return { sol::UNBOUNDED_NO_FEAS, "unbounded, no feasible solution returned" };
  case 4:
    if (timFileData_.barstatus == 1)
      return { sol::UNCERTAIN, "numerical difficulties but possibly optimal" };
  default:
    ;
  }
  switch (timFileData_.barstatus)
  {
  case 2:
    return { sol::LIMIT_NO_FEAS_NODES, "node limit, without a feasible solution" };
  case 3:
    return { sol::LIMIT_NO_FEAS_ITER, "iteration limit, without a feasible solution" };
  case 4:
    return { sol::LIMIT_NO_FEAS_TIME, "time limit, without a feasible solution" };
  case 5:
    return { sol::NUMERIC, "terminated due to unrecoverable numerical issues" };
  case 6:
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted, no feasible solution" };
  case 7:
    return { sol::FAILURE, "too little memory" };
  case 9:
    return { sol::FAILURE, "BARON syntax error (should not happen)" };
  case 10:
    return { sol::FAILURE, "licensing error" };
  case 11:
    return { sol::LIMIT_NO_FEAS_INTERRUPT, "interrupted by Control-C" };
  }
  return { sol::UNKNOWN, fmt::format("Unknown solution status: {}", timFileData_.barstatus) };
}


void BaronmpBackend::FinishOptionParsing() {
  initDirectories(stub(), baronOptions().scratch,
    baronOptions().overwrite);
  copy_common_info_to_other();
  set_verbose_mode(baronOptions().outlev > 0);
  // Write header to baron file
  writeBaron("// BARON {}.{}.{} ({}.{}.{})\n", v_year, v_month, v_day, v_year, v_month, v_day);
}


////////////////////////////// OPTIONS /////////////////////////////////


static const mp::OptionValueInfo iismethod_values[] = {
  { "1", "Fast heuristic", 1},
  { "2", "Deletion filtering algorithm", 2},
  { "3", "Addition filtering algorithm", 3},
  { "4", "Addition-deletion filtering algorithm", 4},
  { "5", "Depth-first search algorithm", 5}
};


static const mp::OptionValueInfo iisorder_values[] = {
  { "-1", "Automatic", -1},
  { "1", "Problem order (as in reformulated model)", 1},
  { "2", "Ascending order by degree", 2},
  { "3", "Descending order by degree", 3},
  { ">=4", "Random order with seed iisorder", 4} 
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

void BaronmpBackend::InitCustomOptions() {

  set_option_header(
      "BARONMP Optimizer Options for AMPL\n"
      "--------------------------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to the "
      "AMPL option ``baronmp_options``. For example::\n"
      "\n"
      "  ampl: option baronmp_options 'mipgap=1e-6';\n");



  #define ADDOPTION(var, n)\
  AddStoredOption(n, baronOptions().description(n), baronOptions().var)

  #define ADDALGOPTION(n) ADDOPTION(n, "alg:"#n " "#n)
  #define ADDTECHOPTION(n) ADDOPTION(n,"tech:"#n " "#n)
  
  char* name;

  ADDALGOPTION(deltaa);
  ADDALGOPTION(deltaterm);
  ADDALGOPTION(deltar);
  ADDALGOPTION(deltat);
  ADDALGOPTION(epsa);
  ADDALGOPTION(epsr);
  ADDALGOPTION(firstfeas);
  ADDALGOPTION(firstloc);
  ADDALGOPTION(iisint);
  ADDALGOPTION(iismethod);
  ADDALGOPTION(iisorder);

  ADDTECHOPTION(barstats);
  ADDTECHOPTION(keepsol);
  ADDTECHOPTION(lpsolver);
  ADDTECHOPTION(cplexlibname);
  ADDTECHOPTION(xpresslibname);
  ADDTECHOPTION(lsolmsg);
  ADDTECHOPTION(lsolver);

  ADDOPTION(maxiter, "lim:iter iterlimit maxiter");
  ADDOPTION(maxtime, "lim:time timelim timelimit maxtime");

  ADDOPTION(nlpsolver, "tech:nlpsolver nlpsol nlpsolver");
  
  ADDTECHOPTION(numsol);
  ADDTECHOPTION(objbound);
  ADDTECHOPTION(objno);
  ADDTECHOPTION(optfile);

  
  ADDTECHOPTION(outlev);
  ADDTECHOPTION(overwrite); // TODORemove
  ADDTECHOPTION(prfreq);
  ADDTECHOPTION(prloc);
  ADDTECHOPTION(problem);
  ADDTECHOPTION(prtime);
  ADDTECHOPTION(scratch);
  ADDTECHOPTION(seed);
  ADDTECHOPTION(sumfile);
  ADDTECHOPTION(threads);
  ADDTECHOPTION(trace);

  AddListOption("tech:optionnative optionnative optnative tech:param",
    "General way to specify values of both documented and "
    "undocumented Baron parameters; value should be a quoted "
    "string (delimited by ' or \") containing a parameter name, a "
    "space, and the value to be assigned to the parameter.",
    baronOptions().inlineparams);

  ////////////////// CUSTOM RESULT CODES ///////////////////
  // todo
  AddSolveResults( {
                     { sol::FAILURE+1, "fatal error 1" },
                     { sol::FAILURE+2, "fatal error 2" },
                     { sol::LIMIT_FEAS_NEW + 1, "AI iteration limit, feasible solution" },
                     { sol::LIMIT_NO_FEAS_NEW + 1, "AI iteration limit, no feasible solution" }
                   } );     // No replacement, make sure they are new
}


void BaronmpBackend::ComputeIIS() {
  // TODO
  //BARONMP_CCALL(BARONMP_ComputeIIS(lp()));
  SetStatus(GetSolveResult());   // could be new information
}

IIS BaronmpBackend::GetIIS() {
  auto variis = VarsIIS();
  auto coniis = ConsIIS();
  auto mv = GetValuePresolver().PostsolveIIS(
    { variis, coniis });
  return { mv.GetVarValues()(), mv.GetConValues()() };
}

ArrayRef<int> BaronmpBackend::VarsIIS() {
  return ArrayRef<int>();
}
pre::ValueMapInt BaronmpBackend::ConsIIS() {

  return { {{ 0, std::vector<int>()}} };
}

} // namespace mp


// AMPLs
void* AMPLSOpenBaronmp(
  const char* slv_opt, CCallbacks cb = {}) {
  return AMPLS__internal__Open(
        std::unique_ptr<mp::BasicBackend>{new mp::BaronmpBackend()},
        cb);
}

void AMPLSCloseBaronmp(AMPLS_MP_Solver* slv) {
  AMPLS__internal__Close(slv);
}

void* GetBaronmpmodel(AMPLS_MP_Solver* slv) {
  return
    dynamic_cast<mp::BaronmpBackend*>(AMPLSGetBackend(slv))->lp();
}
