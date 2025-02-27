#ifndef MP_FLAT_SOL_CHECK_H
#define MP_FLAT_SOL_CHECK_H

#include <cmath>
#include <map>

#include "mp/flat/constr_std.h"
#include "mp/flat/constr_keeper.h"
#include "mp/common.h"
#include "mp/valcvt.h"

namespace mp {

/// A mix-in base class
/// for solution checking.
template <class Impl>
class SolutionChecker {
public:
  /// Check postsolved & unpostsolved (with aux vars) solutions
  /// in various ways.
  /// @param p_extra: !=0 means known infeas solution
  bool CheckSolution(
        ArrayRef<double> x,
        const pre::ValueMapDbl& duals,   // unused
        ArrayRef<double> obj,
        void* p_extra) {
    bool fKnownInfeas = (bool)p_extra;
    if (fKnownInfeas && !MPCD( sol_check_infeas() ))
      return true;
    std::string msgreal, msgidea;
    std::string err_msg;
    try {                            // protect
      std::vector<double> x_back = x;
      if (MPCD( sol_check_mode() ) & (1+2+4+8+16)) {
        msgreal = DoCheckSol(x, duals, obj, {}, x_back, false);
      }
      if (MPCD( sol_check_mode() ) & (32+64+128+256+512)) {
        auto x1 = RecomputeAuxVars(x);
        msgidea = DoCheckSol(x1, duals, obj, x_back, x_back, true);
      }
      if (msgreal.size() || msgidea.size()) {     // Checks failed
        std::string warn
            = "Type                         MaxAbs [Name]   MaxRel [Name]\n";
        warn += msgidea;
        if (msgreal.size()) {
          for (auto i=msgreal.size()-1; i--; ) {
            if ('\n'==msgreal[i] && ' '==msgreal[i+1])
              msgreal[i+1] = '*';
          }
          if (' '==msgreal.front())
            msgreal.front() = '*';
          warn += msgreal;
          warn += "*: Using the solver's aux variable values.\n";
        }
        warn += "Documentation: mp.ampl.com/modeling-tools.html#automatic-solution-check.";
        // Should messages for realistic and idealistic
        // modes be more coordinated?
        // I.e., when no expressions.
        // What if this is an intermediate solution?
        // Should be fine - warning by default,
        // fail if requested explicitly.
        // If warning, we should add the report
        // to that solution's solve message, and
        // a summary in the final solve message.
        // For now, do this via warnings?
        if (MPCD( sol_check_fail() ))
          MP_RAISE_WITH_CODE(int(sol::MP_SOLUTION_CHECK),
                             "Solution check failed - reporting as fatal:\n"
                             + warn);
        else
          MPD( AddWarning(
                 MPD( GetEnv() ).GetSolCheckWarningKey(true),
                 warn,
                 true) );  // replace for multiple solutions
      }
    } catch (const mp::Error& err) {
      if (MPCD( sol_check_fail() ))
        throw;
      err_msg = err.what();
    } catch (const std::exception& exc) {
      err_msg = exc.what();
    } catch (...) {
      err_msg = "unknown error";
    }
    if (err_msg.size()) {
      err_msg += '\n' + msgreal + msgidea;
      if (MPCD( sol_check_fail() ))
        MP_RAISE("Solution check aborted: " + err_msg);
      MPD( AddWarning("Solution check aborted", err_msg) );
      return false;
    }
    return msgreal.empty() && msgidea.empty();
  }

private:
  /// Functor to recompute auxiliary var \a i
  VarsRecomputeFn recomp_fn
  = [this](int i, const VarInfoRecomp& x) {
    if (MPCD( HasInitExpression(i) )) {
      const auto& iexpr = MPCD( GetInitExpression(i) );
      auto pCK = iexpr.GetCK();
#if __GNUC__ > 8            // manylinux2010 fails
      auto resvar = pCK->GetResultVar(iexpr.GetIndex());
      assert(resvar == i);
#endif
      if (!pCK->IsUnused(iexpr.GetIndex()))
        return pCK->ComputeValue(iexpr.GetIndex(), x);
    }
    return x.get_x().get_x()[i];  // no recomputation
  };

public:
  /// Recompute auxiliary variables
  /// @param is_final: if provided, which vars are final
  ArrayRef<double> RecomputeAuxVars(
      ArrayRef<double> x, std::vector<bool> is_final={}) {
    VarInfoRecomp vir {
      MPCD( sol_feas_tol() ),
          true,        // currently not relevant for recomputation
      {x, recomp_fn, is_final},
      {},              // no raw values
      MPCD( GetModel() ).var_type_vec(),
          MPCD( GetModel() ).var_lb_vec(),
          MPCD( GetModel() ).var_ub_vec(),
          MPCD( sol_round() ), MPCD( sol_prec() )
    };
    vir.get_x().set_p_var_info(&vir);
    for (auto i=vir.size(); i--; )
      vir[i];         // touch the variable to be recomputed
    return std::move(vir.get_x().get_x());
  }

protected:
  /// Check single unpostsolved solution.
  /// @param x_raw: for idealistic check,
  ///   should contain auxiliary values from the solver.
  /// @param x_back: solution vector from realistic mode.
  /// It can be changed by the :round and :prec options.
  /// Its auxiliary vars are compared
  /// with recomputed expression values.
  std::string DoCheckSol(
      ArrayRef<double> x,
      const pre::ValueMapDbl& duals,
      ArrayRef<double> obj,
      ArrayRef<double> x_raw,
      std::vector<double>& x_back,
      bool if_recomp_vals) {
    SolCheck chk(x, duals, obj, x_raw,
                 MPCD( GetModel() ).var_type_vec(),
                 MPCD( GetModel() ).var_lb_vec(),
                 MPCD( GetModel() ).var_ub_vec(),
                 MPCD( sol_feas_tol() ),
                 MPCD( sol_feas_tol_rel() ),
                 MPCD( sol_round() ),
                 MPCD( sol_prec() ),
                 if_recomp_vals,
                 if_recomp_vals
                 ? (MPCD( sol_check_mode() ) >> 5)
                 : (MPCD( sol_check_mode() )));  // for raw values
    if (chk.check_mode() & 1)
      CheckVars(chk);
    if (chk.check_mode() & (2+4+8))
      CheckCons(chk);
    if (chk.check_mode() & 16)
      CheckObjs(chk);
    MPD( GenerateViolationsReport(chk, if_recomp_vals) );
    x_back = std::move(chk.x_ext().get_x()); // to reuse 'realistic' vector
    return chk.GetReport();
  }

  void CheckVars(SolCheck& chk) {
    for (auto i=MPCD( num_vars() ); i--; ) {
      auto x = chk.x(i);
      bool aux = !MPCD( is_var_original(i) );
      // no aux vars for idealistic mode
      if (!aux || !chk.if_recomputed()) {
        chk.VarViolBnds().at(aux).CheckViol(
              {MPCD( lb(i) ) - x, MPCD( lb(i) )},
              MPCD( sol_feas_tol() ), MPCD( sol_feas_tol_rel() ),
              MPCD( GetModel() ).var_name(i));
        chk.VarViolBnds().at(aux).CheckViol(
              {x - MPCD( ub(i) ), MPCD( ub(i) )},
              MPCD( sol_feas_tol() ), MPCD( sol_feas_tol_rel() ),
              MPCD( GetModel() ).var_name(i));
        if (MPCD( is_var_integer(i) ))
          chk.VarViolIntty().at(aux).CheckViol(
                { std::fabs(x - std::round(x)),
                  std::round(x) },
                MPCD( sol_int_tol() ), INFINITY,
                MPCD( GetModel() ).var_name(i));
      }
    }
  }

  /// Includes logical constraints.
  void CheckCons(SolCheck& chk) {
    MPD( GetModel() ).ComputeViolations(chk);
  }

  void CheckObjs(SolCheck& chk) {
    const auto& objs = MPCD( GetModel() ).get_objectives();
    // Solvers might have dummy obj.
    // Unbounded problems might have no obj value.
    for (auto i
         =std::min(objs.size(), chk.obj_vals().size());
         i--; ) {
      auto val1 = ComputeValue(objs[i], chk.x_ext());
      chk.ObjViols().CheckViol(
            {std::fabs(chk.obj_vals()[i] - val1), val1},
            MPCD( sol_feas_tol() ), MPCD( sol_feas_tol_rel() ),
            objs[i].name());
    }
  }

  void GenerateViolationsReport(
      SolCheck& chk, bool ) {
    fmt::MemoryWriter wrt;
    if (chk.HasAnyConViols()) {
      Gen1Viol(chk.VarViolBnds().at(0), wrt, true,
               "variable bounds");
      Gen1Viol(chk.VarViolBnds().at(1), wrt, true,
               "aux var bounds");
      Gen1Viol(chk.VarViolIntty().at(0), wrt, true,
               "variable integrality");
      Gen1Viol(chk.VarViolIntty().at(1), wrt, true,
               "aux var integrality");
    }
    GenConViol(chk.ConViolAlg(), wrt, 0);
    GenConViol(chk.ConViolLog(), wrt, 1);
    if (chk.HasAnyObjViols()) {
      Gen1Viol(chk.ObjViols(), wrt, true,
               "objective(s)");
    }
    chk.SetReport( wrt.str() );
  }

  /// Generate message about 1 violation.
  /// @param f_max: whether we need to print
  /// the maximal violations.
  void Gen1Viol(
      const ViolSummary& vs, fmt::MemoryWriter& wrt,
      bool f_max, const std::string& type) {
    if (vs.N_) {
      // wrt.write("  {:16} {:<7}", type, vs.N_);
      wrt.write("  {:27}", type);
      auto vmaxabs = Gen1ViolMax(  // true: 1 for logical
            f_max, vs.epsAbsMax_, vs.nameAbs_);
      auto vmaxrel = Gen1ViolMax(
            f_max, vs.epsRelMax_, vs.nameRel_);
      wrt.write("  {:14}", vmaxabs);
      wrt.write("  {:14}", vmaxrel);
      wrt.write("\n");
    }
  }

  /// Stringify 1 maximal violation
  std::string Gen1ViolMax(
      bool f_max, double viol, const char* nm) {
    fmt::MemoryWriter wrt;
    if (viol>0.0) {
      if (f_max)
        wrt.write("{:.0E}", viol);
      if (nm && *nm != '\0') {
        wrt.write(f_max ? " [" : "[");
        wrt.write("{}]", nm);
      }
    }
    if (0 == wrt.size())
      wrt.write("-");
    return wrt.str();
  }

  /// std::string classnm = alg_log ? "Logical" : "Algebraic";
  /// wrt.write(classnm + " expression violations:\n");
  void GenConViol(
      const std::map< std::string, ViolSummArray<3> >& cvmap,
      fmt::MemoryWriter& wrt, int alg_log) {
    if (cvmap.size()) {
      for (const auto& cva: cvmap) {
        Gen1Viol(cva.second.at(0), wrt, !alg_log,
                 0==cva.first.compare(0, 4, "_lin")
                 ? "algebraic con(s)"
                 : 0==cva.first.compare(0, 5, "_quad")
                   ? "quadratic con(s)"
                 : "expr '"
                 + std::string(cva.first) + "'");
        Gen1Viol(cva.second.at(1), wrt, !alg_log,
                 "interm expr '"
                 + std::string(cva.first) + "'");
        Gen1Viol(cva.second.at(2), wrt, !alg_log,
                 "final expr '"
                 + std::string(cva.first) + "'");
      }
    }
  }
};

}  // namespace mp

#endif // MP_FLAT_SOL_CHECK_H
