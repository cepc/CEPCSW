from valprod.workflow.TestWrapper import TestWrapper
import copy

class Workflow:
  '''
    A workflow is a set of test cases that run in a certain sequence.
    The creation of the TestWrapper objects will be delayed to the moment before the test cases are executed.
  '''

  def __init__(self):
    self.overallCFG = None
    self.stepCFG = []
    self.steps = []

  def addStep(self, name, cmd, **kwa):
    # Creation of the TestWrapper objects will be delayed
    self.stepCFG.append([name, cmd.split(), kwa]) 

  def setOverallCFG(self, cfg):
    self.overallCFG = cfg

  def run(self):
    for sc in self.stepCFG:
      cfg = copy.deepcopy(self.overallCFG)
      for k,v in list(sc[2].items()):
        cfg.setAttr(k,v)
      self.steps.append(TestWrapper(sc[0], sc[1], cfg))
    for step in self.steps:
      ok, what = step.run()
      if not ok:
        return False, what
    return True, ''
