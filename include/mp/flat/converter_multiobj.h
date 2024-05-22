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
  bool PrepareMOIteration() {
    switch (status_) {
    case MOManagerStatus::NOT_SET:
      MP_RAISE("FlatConverter: MultiobjManager not started");
    case MOManagerStatus::NOT_ACTIVE:
      MP_RAISE("FlatConverter: MultiobjManager not running");
    case MOManagerStatus::RUNNING:
      return DoPrepareNextMultiobjSolve();
    case MOManagerStatus::FINISHED:
      return false;
    }
    return false;
  }

  /// Process an unpostsolved solution of the current iteration.
  /// Necessary to be called.
  /// @note Can be called before or after the postsolved solution.
  void ProcessMOIterationUnpostsolvedSolution(pre::ModelValuesDbl& sol) {
    if (IsMOActive()) {
      auto& objs = sol.GetObjValues()();
      assert(1 == objs.size());
      objval_last_ = objs.front();           // 0. save emulated obj value
      // @todo 1. check if the solver correctly reports the current emulated obj
      // 2. Let's recompute the original objectives
      objs.resize( MPCD(num_objs()) );
      for (int i=0; i<(int)objs.size(); ++i)
        objs[i] = ComputeValue(MPCD(get_obj(i)), sol.GetVarValues()());
    }
  }

  /// Process a postsolved solution of the current iteration.
  /// Necessary to be called.
  /// @note Can be called before or after unpostsolved solution.
  /// Should contain at least valid solve status.
  void ProcessMOIterationPostsolvedSolution(const Solution& sol, int solst) {
    assert(sol::Status::NOT_SET != solst);
    assert(IsMOActive());
    assert(MOManagerStatus::FINISHED != status_);
    if ((sol::IsProblemSolvedOrFeasible((sol::Status)solst)
        // || sol::IsProblemMaybeSolved(solst)    // Use this?
         ) && !sol::IsProblemUnbounded((sol::Status)solst))
    {}  // continue
    else
      status_ = MOManagerStatus::FINISHED;
    // We ignore the solution here - but having it provided
    // guarantees that the postsolve has been run.
  }


protected:
  void SetupMultiobjEmulation() {
    status_ = MOManager::MOManagerStatus::RUNNING;
    MPD(set_skip_pushing_objs());  // could have a cleaner system of linking
      // via custom link restoring original objective values,
      // instead of current manual postsolving in ValuePresolver::PostsolveSolution().
    obj_new_ = MPD( get_objectives() );   // no linking
    if (MPD( GetEnv() ).verbose_mode())
      MPD( GetEnv() ).Print(
          "\n\n"
          "==============================================================================\n"
          "MULTI-OBJECTIVE MODE: starting with {} objectives ...\n"
          "==============================================================================\n"
          "==============================================================================\n\n"
          , obj_new_.size());
  }

  /// Do prepare next solve
  bool DoPrepareNextMultiobjSolve() {
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
    MPD( GetModelAPI() ).InitProblemModificationPhase(
        MPD( GetModelInfo() ));
    ReplaceCurrentObj();
    if (i_current_obj_)
      RestrictLastObjVal();
    MPD( GetModelAPI() ).FinishProblemModificationPhase();
    return true;
  }

  void RestrictLastObjVal() {
    assert(i_current_obj_ && i_current_obj_<obj_new_.size());
    const auto& obj_last = obj_new_[i_current_obj_-1];
    if (obj_last.GetQPTerms().size()) {
      if (obj::MAX == obj_last.obj_sense())
        MPD( GetModelAPI() ).AddConstraint(
            QuadConGE{ { obj_last.GetLinTerms(), obj_last.GetQPTerms() }, objval_last_ } );
      else
        MPD( GetModelAPI() ).AddConstraint(
            QuadConLE{ { obj_last.GetLinTerms(), obj_last.GetQPTerms() }, objval_last_ } );
    } else {
      if (obj::MAX == obj_last.obj_sense())
        MPD( GetModelAPI() ).AddConstraint(
            LinConGE{ { obj_last.GetLinTerms() }, objval_last_ } );
      else
        MPD( GetModelAPI() ).AddConstraint(
            LinConLE{ { obj_last.GetLinTerms() }, objval_last_ } );
    }
  }

  void ReplaceCurrentObj() {
    assert(i_current_obj_>=0 && i_current_obj_<obj_new_.size());
    MPD( SetObjectiveTo( MPD(GetModelAPI()), 0, obj_new_[i_current_obj_]) );
  }


private:
  MOManagerStatus status_ {MOManagerStatus::NOT_SET};
  std::vector<QuadraticObjective> obj_new_;     // ranked aggregated objectives
  int i_current_obj_ {-1};
  double objval_last_ {};
};

}  // namespace mp

#endif // CONVERTER_MULTIOBJ_H
