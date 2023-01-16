#!/usr/bin/python
# -*- coding:utf-8 -*-

import subprocess
import time, datetime, os
from valprod.utils.monitors import *
from valprod.workflow.Process import Process,status

class MonitoredProcess(Process):

  def __init__(self, name, cmd, cfg):
    Process.__init__(self, name, cmd, cfg)
    self.pid = None

    # Set up monitoring time interval
    self.interval = self.cfg.getAttr('timeInterval')
    assert type(self.interval) == int or type(self.interval) == float, 'attribute time interval must be a number'
    assert self.interval <= 60 and self.interval >= 0.3, 'attribute time interval must be a number between 0.3 and 60'

    if self.cfg.getAttr('monitorBackend') not in ['ps','prmon']:
      assert False, 'Unknown monitor backend: %s' % self.cfg.getAttr('monitorBackend')

    self.stdout = self.stderr = open(self.logFileName, 'wb+')

  def run(self):
    self._start_process()

    # Create monitors
    all_pids = [self.parent_pid] + self.child_pids
    monitorBackend = self.cfg.getAttr('monitorBackend')
    if monitorBackend == 'ps':
      self._monitor = PSMonitor(self.name, self.interval, all_pids, self.cfg.getAttr('plottingBackend'))
    elif monitorBackend == 'prmon': 
      self._monitor = PrmonMonitor(self.name, self.interval, self.parent_pid)
    if self.cfg.getAttr('perf'): 
      self._perf_monitor = PerfMonitor(self.name, self.parent_pid)
      self._perf_monitor.initialize()
    if self.cfg.getAttr('profileIO'):
      self._io_monitor = DiskIOMonitor(self.name, self.interval, all_pids, self.cfg.getAttr('plottingBackend'))
      self._io_monitor.initialize()
    self._monitor.initialize()

    # The main loop. Wait until the process finishes or the limit is reached
    while True:
      time.sleep(self.interval)
      if not self.process.poll() == None:
        break
      self._monitor.execute()
      if self.cfg.getAttr('profileIO'):
        self._io_monitor.execute()
      if not self._checkLimit():
        break

    self._monitor.finalize()
    if self.cfg.getAttr('perf'):     
      self._perf_monitor.finalize()
    if self.cfg.getAttr('profileIO'):
      self._io_monitor.finalize()
    
    ## Check if the process ends successfully
    self._burnProcess()
    if self.status == status.SUCCESS and self.name:
      self.stdout.close()
      self._parseLogFile()
    if not self.genLog:
      os.remove(self.logFileName)

  def _parseLogFile(self):
    if self.logParser:
      result, self.fatalLine = self.logParser.parseFile(self.logFileName)
      if not result:
        self.status = status.FAIL
