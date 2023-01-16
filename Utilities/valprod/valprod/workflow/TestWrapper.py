import copy
import sys
from valprod.workflow.MonitoredProcess import MonitoredProcess
from valprod.workflow.Process import Process
from valprod.utils.TestConfig import *

class TestWrapper:

  def __init__(self, name, cmd, cfg=None, **kwa):
    self.subprocess = None
    if cfg:
      self.cfg = cfg
    else:
      self.cfg = copy.deepcopy(globalConfig)
    self.cfg.update(**kwa)

    if self.cfg.getAttr('profile') and (not self.cfg.getAttr('monitorBackend')):
      # Ues ps by default
      self.cfg.setAttr('monitorBackend', 'ps')
      

    # Setup sub-process
    if self.cfg.getAttr('monitorBackend'):
      self.subprocess = MonitoredProcess(name, cmd, self.cfg)
    else:
      self.subprocess = Process(name, cmd, self.cfg)

    # Setup plot reference
    self.plotTester = None
    plotRef = self.cfg.getAttr('plotRef')
    if plotRef:
      assert type(plotRef) == str or type(plotRef) == list
      from valprod.utils.PlotTester import PlotTester
      self.plotTester = PlotTester(self.cfg, plotRef)

  def run(self):
    self.subprocess.run()
    ok, summary = self.subprocess.outcome()
    # If process ends succefully, invoke plotTester, if there's one
    if ok and self.plotTester:
      ok, dec = self.plotTester.run()
    return ok, summary
