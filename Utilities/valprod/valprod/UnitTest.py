'''
  example: 
       mytest = UnitTest()
       mytest.addCase('test1', 'python run.py')
       mytest.run()

  1. overall options:

    setTimeLimit: 
      Time limit (s) of the test case. If running time is exceeded, test case will be killed. (default: None)
    setVIRtLimit:
      Virtual memory limit (mb) of the test case. If exceeded, test case will be killed. (default: None)
    enableCPUMonitor/enableVIRMonitor/enableRESMonitor:
      If CPU monitor/ virtual memory monitor/ resident memory monitor is enabled, the real time CPU rate/ virtual memory/ resident memory of the test case will be recorded.
    setFatalPattern:
      If any of the re pattern is detected, in the stdout and stderr of the test case, it will be killed.
      default: {}

  2. local options:
    options of one single test case, example: mytest.addCase('test2', 'python run.py', timeLimit=900, genLog=True, RESMonitor=True)

    genLog:
      If genLog is set to True, the stdout and stderr of the test case will be recorded. (default: False)
    maxTime:
      Time limit (s) of the test case. If running time is exceeded, test case will be killed. (default: None)
    maxVIR:
      Virtual memory limit (kb) of the test case. If exceeded, test case will be killed. (default: None)
    CPUMonitor/RESMonitor/VIRMonitor:
      If CPU monitor/ virtual memory monitor/ resident memory monitor is enabled, the real time CPU rate/ virtual memory/ resident memory of the test case will be recorded. (default: False)
    timeInterval:
      Monitoring interval (s) of the test case. (default: 0.5, range: 0.3-10)
    plotRef:
      The plot reference file. If set, the plot output of the test case will be compared. (default: None)
    histTestMeth:
      Method of the histogram testing. (default: Kolmogorov, Kolmogorov or Chi2)
    histTestCut:
      The P-Value cut of the histogram test. If result is smaller than this cut, the histogram will be treated as inconsistent. (default: 0.9)

  3. workflow:
'''

import unittest
from valprod.workflow.TestWrapper import TestWrapper
from valprod.utils.TestConfig import TestConfig

class UnitTest:

  def __init__(self):
    self.defaultCFG = TestConfig()
    self.case = ValprodTestCase
    self.suite = unittest.TestSuite()
    self.runner = unittest.TextTestRunner()
    self.workflowID = 0

  def addCase(self, name, cmd, **kwa):
    setattr(self.case, "test_" + name, genTestCase(name, cmd, self.defaultCFG, **kwa))
    self.suite.addTest(self.case("test_" + name))

  def addWorkflow(self, workflow):
    self.workflowID += 1
    wfName = "test_workflow%d" % self.workflowID
    setattr(self.case, wfName, genWorkflow(workflow, self.defaultCFG))
    self.suite.addTest(self.case(wfName))

  def run(self):
    self.runner.run(self.suite)

  def setTimeLimit(self, value):
    self.defaultCFG.setAttr("maxTime", value)

  def setVIRLimit(self, value):
    self.defaultCFG.setAttr("maxVIR", value)

  def enableProfile(self):
    self.defaultCFG.setAttr("profile", True)

  def enablePerf(self):
    self.defaultCFG.setAttr("perf", True)

  def enableIOProfile(self):
    self.defaultCFG.setAttr("profileIO", True)

  def setFatalPattern(self, value):
    self.defaultCFG.setAttr("fatalPattern", value)

  def enableGenLog(self):
    self.defaultCFG.setAttr("genLog", True)


class ValprodTestCase(unittest.TestCase):

  def setUp(self):
    pass

  def tearDown(self):
    pass

def genTestCase(name, cmd, cfg, **kwa):
  def case(self):
    run = TestWrapper(name, cmd, cfg, **kwa)
    ok, what = run.run()
    self.assertTrue(ok, what)
  return case

def genWorkflow(workflow, cfg):
  def wkCase(self):
    workflow.setOverallCFG(cfg)
    ok, what = workflow.run()
    self.assertTrue(ok, what)
  return wkCase
