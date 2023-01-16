import datetime
import os,sys
import select
import subprocess
from valprod.utils.shellUtil import GetMemUse
from valprod.utils.Parser import Parser

class status:
  (SUCCESS, FAIL, TIMEOUT, OVERFLOW, ANR) = range(0, 5)
  Description = {
          FAIL: 'Return code is not zero',
          TIMEOUT: 'Run time exceeded',
          OVERFLOW: 'Memory overflow',
          ANR: 'Not responding'
        }
  @staticmethod
  def describe(stat):
    if stat in status.Description:
      return status.Description[stat]
    

class Process:

  def __init__(self, name, exe, cfg):

    print('')
    self.cfg = cfg
    self.name = name or 'test'
    self.logParser = cfg.getAttr('parser') and Parser(cfg) or None
    self.executable = exe
    self.genLog = self.cfg.getAttr('genLog')

    # Log file name depends on what we are running
    logname = self.cfg.getAttr('logName')
    if logname:
      self.logFileName = logname
    else:
      self.logFileName = self.name + '.log'
    if self.cfg.getAttr('step'):
        self.logFileName = self.cfg.getAttr('step')+'.log'

    # Merge stdout and stderr
    self.stdout = self.stderr = subprocess.PIPE
    self.process = None
    self.parent_pid = None
    self.child_pids = []
    self.returncode = None
    self.status = None
    self.memory_peak = 0
    self.timeout = self.cfg.getAttr('timeout')
    self.timeLimit = self.cfg.getAttr('maxTime')
    self.memoryLimit = self.cfg.getAttr('maxMEM')
    self.duration = None
    self.start = None
    self.killed = None
    self.fatalLine = None

  def _start_process(self):
    print('Running test: %s' % self.name)
    print(self.executable)
    self.start = datetime.datetime.now()
    self.process = subprocess.Popen(args = self.executable, stdout = self.stdout, stderr = subprocess.STDOUT)
    self.parent_pid = self.process.pid

    ## Get children pids recursively. This is important in case the parent process
    ## is just a thin wrapper.
    try:
      import psutil
    except:
      return
    try:
      parent = psutil.Process(self.parent_pid)
      children = parent.children(recursive=True)
      for child in children:
        self.child_pids.append(child.pid)
    except psutil.NoSuchProcess:
      ## Parent process is not created or has already exited. Do nothing here
      pass

  def run(self):
    self._start_process()
    if self.genLog:
      # TODO
      # * allow specify log directory
      # * should support directory structure

      logDir = os.path.dirname(self.logFileName)
      if len(logDir) and not os.path.exists(logDir):
        os.makedirs(logDir)
      logFile = open(self.logFileName, 'w')
    stdout_wait_time = 0
    while True:
      stdout_start = datetime.datetime.now()
      fs = select.select([self.process.stdout], [], [], 10)
      stdout_end = datetime.datetime.now()
      if not fs[0]:
        # No response
        stdout_wait_time = stdout_wait_time + (stdout_end - stdout_start).seconds
        if self.timeout and stdout_wait_time > self.timeout:
          self.status = status.ANR
          self._kill()
          break
      else:
        stdout_wait_time = 0 
      if self.process.stdout in fs[0]:
        # Incoming message to parse
        data = os.read(self.process.stdout.fileno(), 1024).decode()
        if not data:
          break
        # If it is called in analysis step, we print the log info to screen
        if self.cfg.getAttr("step"):
          for l in data.splitlines(): print("[%d]: "%self.parent_pid, l)
        if self.genLog:
          logFile.write(str(data)+"\n\n")
        if self.logParser:
          if not self._parseLog(data):
            self.status = status.FAIL
            self._kill()
            break
      if not self._checkLimit():
        break
    self._burnProcess()
    self._checkLimit()
    if self.genLog:
      logFile.close()

  def getDuration(self):
    return self.duration

  def _checkLimit(self):
    self.duration = (datetime.datetime.now() - self.start).seconds
    if self.timeLimit and (self.duration >= self.timeLimit):
      # Time out
      self.status = status.TIMEOUT
      self._kill()
      return False
    if self.memoryLimit and (self._getMem() >= self.memoryLimit):
      # Memory overflow
      self.status = status.OVERFLOW
      self._kill()
      return False
    return True

  def _kill(self):
    if not self.process:
      return
    import signal
    try:
      os.kill(self.parent_pid, signal.SIGKILL)
      os.waitpid(-1, os.WNOHANG)
    except:
      pass
    for child_pid in self.child_pids:
      try:
        os.kill(child_pid, signal.SIGKILL)
      except:
        pass

  def _parseLog(self, data):
    result, self.fatalLine = self.logParser.parse(data)
    return result

  def _burnProcess(self):
    self.returncode = self.process.wait()
    if self.status:
      return
    if 0 == self.returncode:
      self.status = status.SUCCESS
    else:
      self.status = status.FAIL
    #FIXME: it seems that root macro process won't give a 0 return code
    if type(self.executable) == list and self.executable[0] == 'root':
      self.status = status.SUCCESS
    if type(self.executable) == str and self.executable.startswith('root'):
      self.status = status.SUCCESS

  def _getMem(self):
    if not self.parent_pid:
      return 0
    else:
      mem_now = GetMemUse(self.parent_pid)
      for child_pid in self.child_pids:
        mem_now = mem_now + GetMemUse(child_pid)
      if mem_now > self.memory_peak:
        self.memory_peak = mem_now
      return mem_now

  def outcome(self):
    summary = ''
    if self.duration:
      summary = 'Running time: %ds. ' % self.duration 
    if self.memory_peak:
      summary = summary + 'Peak memory: %dMb' % self.memory_peak
    if self.status == status.SUCCESS:
      summary = 'Successful. \n' + summary
      return True, summary
    if self.fatalLine:
      summary = 'Failed. \n' + 'FatalLine: ' + self.fatalLine + '\n' + summary
      return False, summary
    summary = 'Failed. \n' + status.describe(self.status) + '\n' + summary
    return False, summary
