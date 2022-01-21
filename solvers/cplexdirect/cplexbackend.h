#ifndef MP_CPLEX_BACKEND_H_
#define MP_CPLEX_BACKEND_H_

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#elif _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

extern "C" {
  #include <ilcplex/cplex.h>
}

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include <string>

#include "mp/flat/backend.h"
#include "mp/flat/std_constr.h"

namespace mp {

class CplexBackend :
    public Backend<CplexBackend>
{
  using BaseBackend = Backend<CplexBackend>;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  CplexBackend();
  ~CplexBackend();

  /// Metadata
  static const char* GetSolverName() { return "IBM ILOG CPLEX"; }
  std::string GetSolverVersion();
  static const char* GetSolverInvocationName();
  static const char* GetSolverLongName() { return nullptr; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  /// [[ Prototype the incremental interface ]]
  void InitProblemModificationPhase();
  void FinishProblemModificationPhase();

  void AddVariables(const VarArrayDef& );
  void SetLinearObjective( int iobj, const LinearObjective& lo );

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseBackend)

  ACCEPT_CONSTRAINT(RangeLinCon, Recommended, CG_Linear)
  void AddConstraint(const RangeLinCon& lc);
  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, AcceptedButNotRecommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);


  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const;
  bool IsQCP() const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;

  void ExportModel(const std::string& file);


  //////////////////////////// SOLVING ///////////////////////////////
  void SetInterrupter(mp::Interrupter* inter);
  void SolveAndReportIntermediateResults();

  Solution GetSolution();
  double ObjectiveValue() const;

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism
  void OpenSolver();
  void CloseSolver();
  void InitCustomOptions();

  static double Infinity() { return CPX_INFBOUND; }
  static double MinusInfinity() { return -CPX_INFBOUND; }

protected:
  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution();
  pre::ValueMapDbl DualSolution();
  ArrayRef<double> DualSolution_LP();

  void WindupCPLEXSolve();

  std::pair<int, std::string> ConvertCPLEXStatus();
  void AddCPLEXMessages();

private:
  CPXENVptr     env = NULL;
  CPXLPptr      lp = NULL;

  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
  };
  Options storedOptions_;

public:
  /// These methods access CPLEX options. Used by AddSolverOption()
  void GetSolverOption(int key, int& value) const;
  void SetSolverOption(int key, int value);
  void GetSolverOption(int key, double& value) const;
  void SetSolverOption(int key, double value);
  void GetSolverOption(int key, std::string& value) const;
  void SetSolverOption(int key, const std::string& value);

};

}  // namespace mp

#endif  // MP_CPLEX_BACKEND_H_
