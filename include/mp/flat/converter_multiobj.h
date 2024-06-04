#ifndef CONVERTER_MULTIOBJ_H
#define CONVERTER_MULTIOBJ_H

#include <vector>
#include <cassert>

#include "mp/env.h"
#include "mp/common.h"
#include "mp/valcvt-base.h"
#include "mp/error.h"
#include "mp/flat/obj_std.h"

namespace mp {

/// A mix-in base class managing multiobjective emulation
template <class Impl>
class MOManager {
protected:
  /// MOManager status
  enum class MOManagerStatus {
    NOT_SET,
    NOT_ACTIVE,
    RUNNING,
    FINISHED
  };

  /// See if we need to emulate multiple objectives.
  /// @note This should be called first.
  void ConsiderEmulatingMultiobj() {
    status_ = MOManagerStatus::NOT_ACTIVE;
    if (MPCD(num_objs())>1       // have multiple objectives
        && MPCD(GetEnv()).multiobj_has_native()==false)
      SetupMultiobjEmulation();
    // Anything todo otherwise?
  }


public:
  /// Is MOManager active?
  /// This is relevant after initialization via
  /// ConsiderEmulatingMultiobj().
  /// @note Check this before using other public methods.
  bool IsMOActive() const {
    return
        MOManagerStatus::RUNNING==status_
           || MOManagerStatus::FINISHED==status_;
  }

  /// Prepare next multiobj iteration?
  /// @note Call this before a MO iteration.
  /// @return true iff the model is ready for the next iteration.
  bool PrepareMOIteration(
      std::function<sol::Status(void)> get_stt, std::function<Solution(void)> get_sol) {
    switch (status_) {
    case MOManagerStatus::NOT_SET:
      MP_RAISE("FlatConverter: MultiobjManager not set up");
    case MOManagerStatus::NOT_ACTIVE:
      MP_RAISE("FlatConverter: MultiobjManager not running");
    case MOManagerStatus::RUNNING:
      return DoPrepareNextMultiobjSolve(get_stt, get_sol);
    case MOManagerStatus::FINISHED:
      return false;
    }
    return false;
  }

  /// Obtain and process a postsolved solution of the current iteration.
  /// Can implicitly call ProcessMOIterationUnpostsolvedSolution().
  /// @return (whether we should continue, solve_result).
  std::pair<bool, sol::Status> ProcessMOIterationPostsolvedSolution(
      std::function<sol::Status(void)> get_stt, std::function<Solution(void)> get_sol) {
    auto solst = get_stt();
    assert(sol::Status::NOT_SET != solst);
    assert(IsMOActive());
    assert(MOManagerStatus::FINISHED != status_);
    if ((sol::IsProblemSolvedOrFeasible((sol::Status)solst)
        // || sol::IsProblemMaybeSolved(solst)    // Use this?
         )  && !sol::IsProblemUnbounded((sol::Status)solst)  // Don't want unbounded
        ) {                                                  // (but LIMIT can have this undiscovered)
      get_sol();    // This implicitly calls ProcessMOIterationUnpostsolvedSolution().
      return { true, solst };
      // We ignore the solution here - but having it provided
      // guarantees that the postsolve has been run.
    }
    status_ = MOManagerStatus::FINISHED;
    return { false, solst };
  }

  /// Process an unpostsolved solution of the current iteration.
  /// This is called from solution postsolve initiated by solution getter.
  /// @note Can be called before or after the postsolved solution.
  void ProcessMOIterationUnpostsolvedSolution(pre::ModelValuesDbl& sol) {
    if (IsMOActive()) {
      auto& objs = sol.GetObjValues()();
      if (objs.size())
        objval_last_ = objs.front();           // 0. save emulated obj value
      // @todo 1. check if the solver correctly reports the current emulated obj
      // 2. Let's recompute the original objectives
      const auto& xx = sol.GetVarValues()();
      if (xx.size()) {                         // This can be invoked w/o solution
        objs.resize( MPCD(num_objs()) );
        for (int i=0; i<(int)objs.size(); ++i)
          objs[i] = ComputeValue(MPCD(get_obj(i)), xx);
      }
    }
  }


protected:
  void SetupMultiobjEmulation() {
    status_ = MOManager::MOManagerStatus::RUNNING;
    MPD(set_skip_pushing_objs());  // could have a cleaner system of linking
      // via custom link restoring original objective values,
      // instead of current manual postsolving in ValuePresolver::PostsolveSolution().
    const auto& obj_orig = MPD( get_objectives() );   // no linking
    ///////////////// Read / set default suffixes ///////////////////
    std::vector<int> objpr = MPD( ReadIntSuffix( {"objpriority", suf::OBJ} ) );  // int only
    objpr.resize(obj_orig.size(), 0.0);               // blend objectives by default
    std::vector<double> objwgt = MPD( GetMOWeightsLegacy() );
    if (objwgt.empty()) {
      objwgt.resize(obj_orig.size(), 1.0);            // Default "intuitive" weights
      FlipDiffSenseSigns(objwgt);                     // We handle "legacy" format below
    }
    std::vector<double> objtola = MPD( ReadDblSuffix( {"objabstol", suf::OBJ} ) );
    objtola.resize(obj_orig.size(), 0.0);
    std::vector<double> objtolr = MPD( ReadDblSuffix( {"objreltol", suf::OBJ} ) );
    objtolr.resize(obj_orig.size(), 0.0);
    std::map<int, std::vector<int>, std::greater<int> > pr_map;      // Decreasing order
    for (int i=0; i<objpr.size(); ++i)
      pr_map[objpr[i]].push_back(i);
    obj_new_ = {};         ////////////////// Aggregate new objectives ///////////////////
    obj_new_.reserve(pr_map.size());
    obj_new_tola_.reserve(pr_map.size());
    obj_new_tolr_.reserve(pr_map.size());
    for (const auto& pr_level: pr_map) {
      const auto& i0_vec = pr_level.second;
      obj_new_.push_back(obj_orig.at(i0_vec.front()));
      obj_new_.back().set_sense(obj_orig.front().obj_sense());      // "Legacy" obj:multi:weight
      obj_new_.back().GetLinTerms() *= objwgt.at(i0_vec.front());   // Use weight
      obj_new_.back().GetQPTerms() *= objwgt.at(i0_vec.front());
      obj_new_tola_.push_back(objtola.at(i0_vec.front()));
      obj_new_tolr_.push_back(objtolr.at(i0_vec.front()));
      for (auto i0i=i0_vec.size(); --i0i; ) {
        // Add next objective with weight and sense factor
        double sensef
            = (obj_orig.at(i0_vec.front()).obj_sense() == obj_orig.at(i0_vec[i0i]).obj_sense())
                  ? 1.0 : -1.0;
        auto lt1 = obj_orig.at(i0_vec[i0i]).GetLinTerms();
        lt1 *= sensef * objwgt.at(i0_vec[i0i]);
        auto qt1 = obj_orig.at(i0_vec[i0i]).GetQPTerms();
        qt1 *= sensef * objwgt.at(i0_vec[i0i]);
        obj_new_.back().GetLinTerms().add(lt1);
        obj_new_.back().GetQPTerms().add(qt1);
        // Max the tolerances
        obj_new_tola_.back() = std::max(obj_new_tola_.back(), objtola.at(i0_vec[i0i]));
        obj_new_tolr_.back() = std::max(obj_new_tolr_.back(), objtolr.at(i0_vec[i0i]));
      }
      obj_new_.back().GetLinTerms().sort_terms();
      obj_new_.back().GetQPTerms().sort_terms();
    }
    if (MPD( GetEnv() ).verbose_mode())
      MPD( GetEnv() ).Print(
          "\n\n"
          "==============================================================================\n"
          "MULTI-OBJECTIVE MODE: starting with {} objectives ({} combined) ...\n"
          "==============================================================================\n"
          "==============================================================================\n\n"
          , obj_orig.size(), obj_new_.size());
  }

  /// Do prepare next solve
  bool DoPrepareNextMultiobjSolve(
      std::function<sol::Status(void)> get_stt, std::function<Solution(void)> get_sol) {
    if (++i_current_obj_ >= obj_new_.size()) {
      status_ = MOManagerStatus::FINISHED;
      MPD( GetEnv() ).Print(
          "\n\n"
          "==============================================================================\n"
          "MULTI-OBJECTIVE MODE: done.\n\n");
      return false;                              // all done
    }
    if (MPD( GetEnv() ).verbose_mode())
      MPD( GetEnv() ).Print(
          "\n\n"
          "MULTI-OBJECTIVE MODE: objective {} (out of {}) ...\n"
          "==============================================================================\n\n"
          , i_current_obj_+1, obj_new_.size());
    if (i_current_obj_) {
      auto proc_sol = ProcessMOIterationPostsolvedSolution(get_stt, get_sol);
      if (!proc_sol.first) {
        if (MPD( GetEnv() ).verbose_mode())
          MPD( GetEnv() ).Print(
              "   ... ABORTING: previous iteration's solve result: {} (code {}.)\n"
              "==============================================================================\n\n"
              , sol::GetStatusName(proc_sol.second), proc_sol.second);
        return false;
      }
      RestrictLastObjVal();
    }
    MPD( FillConstraintCounters( MPD( GetModelAPI() ), *MPD( GetModelInfoWrt() ) ) );   // @todo a hack.
    MPD( GetModelAPI() ).InitProblemModificationPhase(   // For adding the new constraint. @todo a hack.
        MPD( GetModelInfo() ));                          // Ideally Model would notice changes and notify
    ReplaceCurrentObj();                  // After allowing model modification (needed by SCIP.)
    MPD( AddUnbridgedConstraintsToBackend( MPD( GetModelAPI() ), nullptr) );
    MPD( GetModelAPI() ).FinishProblemModificationPhase();            // ModelAPI automatically.
    return true;
  }

  void RestrictLastObjVal() {
    assert(i_current_obj_ && i_current_obj_<obj_new_.size());
    const auto& obj_last = obj_new_[i_current_obj_-1];
    auto lim = objval_last_;
    auto diff = std::max(                 // Apply degradation tolerance
        obj_new_tola_[i_current_obj_-1], std::fabs(lim) * obj_new_tolr_[i_current_obj_-1]);
    lim += diff * (obj::MAX==obj_last.obj_sense() ? -1.0 : 1.0);
    if (obj_last.GetQPTerms().size()) {
      if (obj::MAX == obj_last.obj_sense())
        MPD( AddConstraint(
            QuadConGE{ { obj_last.GetLinTerms(), obj_last.GetQPTerms() }, objval_last_ } ) );
      else
        MPD( AddConstraint(
            QuadConLE{ { obj_last.GetLinTerms(), obj_last.GetQPTerms() }, objval_last_ } ) );
    } else {
      if (obj::MAX == obj_last.obj_sense())
        MPD( AddConstraint(
            LinConGE{ { obj_last.GetLinTerms() }, objval_last_ } ) );
      else
        MPD( AddConstraint(
            LinConLE{ { obj_last.GetLinTerms() }, objval_last_ } ) );
    }
  }

  void ReplaceCurrentObj() {
    assert(i_current_obj_>=0 && i_current_obj_<obj_new_.size());
    MPD( SetObjectiveTo( MPD(GetModelAPI()), 0, obj_new_[i_current_obj_]) );
  }

  ArrayRef<double> GetMOWeightsLegacy() {
    std::vector<double> objw = MPD( ReadDblSuffix( {"objweight", suf::OBJ} ) );
    if (objw.size() && 2==MPD( GetEnv() ).multiobj_weight()) {   // user wants "intuitive"
      FlipDiffSenseSigns(objw);            // Backend / Emulator want "legacy"
    }
    return objw;
  }

  /// Convert between the options of obj:multi:weight
  void FlipDiffSenseSigns(std::vector<double>& objw) {
    const auto& obj = MPD( get_objectives() );
    for (auto i=obj.size(); --i; )             // forall i>1
      if (obj[i].obj_sense() != obj.front().obj_sense())
        objw[i] = -objw[i];
  }

private:
  MOManagerStatus status_ {MOManagerStatus::NOT_SET};
  std::vector<QuadraticObjective> obj_new_;     // ranked aggregated objectives
  std::vector<double> obj_new_tola_;
  std::vector<double> obj_new_tolr_;
  int i_current_obj_ {-1};
  double objval_last_ {};
};

}  // namespace mp

#endif // CONVERTER_MULTIOBJ_H
