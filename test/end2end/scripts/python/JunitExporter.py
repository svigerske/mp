from junitparser import TestCase, TestSuite, JUnitXml,  Failure, Skipped, Error
from Exporter import Exporter
import Solver
import ModelRunner
from os.path import splitext
import xml.sax.saxutils

class JunitExporter(Exporter):

    def __init__(self, fileName=""):
        super().__init__()
        self._fileName=fileName
        self._test_suites={}
        self._versions={}

    def get_file_name(self, solver: str):
        base_name, _= splitext(self._fileName)
        return f"{solver}-{base_name}.xml"
    
    def _exportInstanceResults(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        if i == 1:
            self.initialize(mr)
        self.append_last_results(mr)
        self.save()
        
    def initialize(self, mr: ModelRunner):
        runs=mr.getRuns()
        runners = mr.getRunners()
        solvers_runs = zip(runners, runs)
        for (runner, solver_run) in solvers_runs:
            last_run=solver_run[-1]
            solver=last_run["solver"]
            if isinstance(solver, Solver.Solver):
                solvername=solver.getName()
            else:
                solvername=solver
                solver=runner
            name = f"{solvername}-{Exporter.get_platform_string()}"
            self._test_suites[solvername]= TestSuite(name)
            
        
    def append_last_results(self, mr: ModelRunner):
      i = len( mr.getRuns()[0] )
      m = mr.getModels()[i-1]
      res = [m.getName(), m.getExpectedObjective()]
      for r in mr.getRuns():
            last_run=r[-1]
            solver=last_run["solver"]
            if isinstance(solver, Solver.Solver):
                solver=solver.getName()

            time= {"solutionTime" : r[-1].get("solutionTime", 0)}
            if "times" in r[-1]:
                for p in ["setup", "solver"]:  
                 time[p] = r[-1]["times"].get(p, 0)

            tc=TestCase(res[0], time=time["solutionTime"])
            if "Skipped" in last_run["outmsg"]:
                tc.result=[Skipped(last_run["outmsg"])]
                res="Skipped"
            elif "script failure" in last_run["outmsg"]:
                tc.result=[Error(last_run["eval_fail_msg"])]
                res="Failure"
            elif "eval_fail_msg" in last_run:
                safe_string = xml.sax.saxutils.escape(last_run["output"])
                safe_string =safe_string.replace('\b', '')
                tc.system_out=safe_string
                tc.result=[Failure(last_run["eval_fail_msg"])]
                res="Failure"
            else:
                res = "OK"
                
            self._test_suites[solver].add_testcase(tc)
            
    def save(self):
      for solver,ts in self._test_suites.items():
        xml = JUnitXml()
        xml.add_testsuite(ts)
        xml.write(self.get_file_name(solver), pretty=True)

            
