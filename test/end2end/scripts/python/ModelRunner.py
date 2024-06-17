from Solver import Solver
from AMPLRunner import AMPLRunner
from TimeMe import TimeMe
from Model import Model

class ModelRunner(object):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, runners, optionsExtra=None):
        self._runners = runners
        self._amplRunners = None
        self._runs = [ list() for r in self._runners ]
        self._optionsExtra = optionsExtra

    def getRuns(self):
        return self._runs

    def getModels(self) -> list:
        return self._models
    def getRunners(self)->list:
        return self._runners
    def getLogFileName(m : Model, s : Solver):
        return f"{m.getName()}.{s.getName()}.log"

    def run(self, modelList: list, exporter=None, keepLogs = False, verbose=False):
        """Run the models in this instance. If exporter != None, it exports the results as it goes"""
        self._models = modelList
        n = 0
        nFailedSolver = [0 for r in self.getRunners()]
        nFailedScriptOrAMPL = [0 for r in self.getRunners()]
        nSkipped = [0 for r in self.getRunners()]
        EFM ="eval_fail_msg"
        EM ="errormsg"
        # Set to true for junit to add the solver output to the test result
        keep_output=False
        try:
        # Attempt to access the class
            __import__("JunitExporter")
            import JunitExporter
            if isinstance(exporter, JunitExporter.JunitExporter):
                keep_output=True
        except:
            pass
        
        for m in modelList:
            n += 1
            if m.isNL():
                cr = self._runners
                msg = "{}. Solving as NL: '{}'".format(n, m.getName())
            else:
                if not self._amplRunners:
                    self._amplRunners = [
                        AMPLRunner(r, self._optionsExtra, 
                                   printOutput=verbose, storeOutput=keep_output)
                        for r in self._runners ]
                cr = self._amplRunners
                msg = "{}. Solving with AMPL: '{}'".format(n, m.getName())
            print("{0: <80}".format(msg), end="", flush=True)
            for (i,r) in enumerate(cr):
                t = TimeMe()
                with t:
                  try:
                    r.isBenchmark = len(cr) > 1
                    self._runs[i].append({})
                    if isinstance(r, AMPLRunner):
                        ss = r.getSolver()
                    else:
                        ss = r
                    if not m.isSubsetOfTags(ss.getSupportedTags()) or \
                            m.hasAnyTag(ss.getUnsupportedTags()):
                        self._runs[i][-1]["outmsg"] = "Skipped, unsupported tags"
                        self._runs[i][-1]["solver"] = ss
                        if r.isBenchmark:
                            print("\n\t\t{0: <20}: Skipped due to unsupported tags".
                              format(ss.getName()), flush=True)
                        else:
                            print("Skipped due to unsupported tags", flush=True)
                        nSkipped[i] += 1
                    else:
                      if keepLogs:
                        r.runAndEvaluate(m, logFile=ModelRunner.getLogFileName(m, ss))
                      else:
                        r.runAndEvaluate(m, logFile=None)
                    
                      stats = r.getSolutionStats()
                      
                      if EFM in stats:
                        if EM in stats:
                            stats[EFM]= stats[EM]
                      if keep_output:
                         stats["output"]=r.get_output()
                      self._runs[i][-1] = stats
                      if exporter:
                        if not exporter.printStatus(m, stats):
                            nFailedSolver[i] += 1
                  except Exception as exc:
                    self._runs[i][-1]["outmsg"] = "AMPL(PY)/script failure"
                    self._runs[i][-1]["solver"] = ss
                    self._runs[i][-1][EFM]=str(exc)
                    print("   EXCEPTION: ", exc)
                    nFailedScriptOrAMPL[i] += 1
                print("  (%.4fs, %d failed solver, %d failed AMPL(PY)/script, %d skipped)" %
                  (t.interval, nFailedSolver[i], nFailedScriptOrAMPL[i], nSkipped[i]),
                  end="", flush=True)
            if exporter:
                self.export(exporter)
            print("")
        exporter.export()

    def export(self, exporter):
        exporter.exportInstanceResults(self)

