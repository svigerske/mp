#ifndef MP_GCG_BACKEND_H_
#define MP_GCG_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "gcgmpcommon.h"

namespace mp {

class GcgBackend :
    public FlatBackend< MIPBackend<GcgBackend> >,
    public GcgCommon
{
  using BaseBackend = FlatBackend< MIPBackend<GcgBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  GcgBackend();
  ~GcgBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "gcg"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-GCG"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "GCG"; }
  /// Version displayed with -v
  std::string GetSolverVersion();
  /// External libraries displayed with -v
  std::string set_external_libs() override;
  
  /// Name for diagnostic messages
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  /// Init custom driver options, such as outlev, writeprob
  void InitCustomOptions() override;
  /// Chance for the Backend to init solver environment, etc.
  void InitOptionParsing() override { }
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;



  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions. 
  USING_STD_FEATURES;

  /**
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, true)

  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, true)
  // TODO If getting/setting a basis is supported, implement the 
  // accessor and the setter below
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

  /**
  * MIP warm start
  **/
  // TODO If MIP warm start is supported, implement the function below
  // to set a non-presolved starting solution
  ALLOW_STD_FEATURE(MIPSTART, false)
  void AddMIPStart(ArrayRef<double> x0,
									 ArrayRef<int> sparsity) override;


 /**
  * Get MIP Gap
  **/
  // TODO Implement to return MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // TODO Implement to return the best dual bound value
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;

  /**
  * Compute the IIS and obtain relevant values
  **/
  ALLOW_STD_FEATURE(IIS, false)
  /// Compute IIS
  void ComputeIIS() override;
  /// Retrieve IIS elements
  IIS GetIIS() override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override;
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

public:  // public for static polymorphism
  /// Solve, no model modification any more (such as feasrelax).
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise/finally via ReportResults()
  void Solve() override;

  /// Default impl of GetObjValues()
  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; } 


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:
  void ExportModel(const std::string& file);

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupGCGSolve();

  void ReportResults() override;
  void ReportGCGResults();

  void ReportGCGPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddGCGMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int>);
  void ConStatii(ArrayRef<int>);

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();


private:
  void InputDecomposition();
  
  /// These options are stored in the class
  struct Options {
    std::string exportFile_, logFile_, paramRead_;
    int heuristics_ = 0;
    int cuts_ = 0;
    int presolvings_ = 0;
    int outlev_ = 0;
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_GCG_BACKEND_H_
