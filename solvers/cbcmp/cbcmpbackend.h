#ifndef MP_CBCMP_BACKEND_H_
#define MP_CBCMP_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cbcmpcommon.h"

namespace mp {




class CbcmpBackend :
    public FlatBackend< MIPBackend<CbcmpBackend> >,
    public CbcmpCommon
{
  using BaseBackend = FlatBackend< MIPBackend<CbcmpBackend> >;
  //////////////////// [[ The public interface ]] //////////////////////
public:
  CbcmpBackend();
  ~CbcmpBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "cbc"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "cbc"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-CBC"; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override { }
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;



  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions. 
  USING_STD_FEATURES;

  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;
  /**
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, false)

  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, true)
  // TODO If getting/setting a basis is supported, implement the 
  // accessor and the setter below
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;


  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE(WARMSTART, true)
  void AddPrimalDualStart(Solution sol0) override;
  /**
  * MIP warm start
  **/
  ALLOW_STD_FEATURE(MIPSTART, true)
	void AddMIPStart(ArrayRef<double> x0,
									 ArrayRef<int> sparsity) override;


 /**
  * Get MIP Gap
  **/
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, false)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // TODO Implement to return the best dual bound value
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, false)
  double BestDualBound() override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override;
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;

  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; } 


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism
  void InitCustomOptions() override;

protected:

  void OpenSolver();
  void CloseSolver();

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupCBCMPSolve();

  void ReportResults() override;
  void ReportCBCMPResults();

  void ReportCBCMPPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddCBCMPMessages();

private:
  /// These options are stored in the class
  struct Options {
    int timeLimit_ = 0;
    int outlev_ = 0;
  };
  Options storedOptions_;
  template <class T> struct OwningOption {
    std::string name;
    std::string desct;
    std::string par;
    T min;
    T max;
  };

  std::vector<OwningOption<int>> int_options;
  std::vector<OwningOption<double>> double_options;

};

}  // namespace mp

#endif  // MP_CBCMP_BACKEND_H_
