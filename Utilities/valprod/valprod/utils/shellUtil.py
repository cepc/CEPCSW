# -*- coding:utf-8 -*-
# author liteng

import time
import os,sys
import subprocess

def GetTotalSwap():
  value = os.popen("free | grep Swap | awk '{print $2}'").read()
  return int(value) / 1024 ## convert from kB to MB

def GetTotalMem():
  value = os.popen("free | grep Mem | awk '{print $2}'").read()
  return int(value) / 1024

def GetTotalCpuRate():
  value = os.popen("vmstat|grep -v procs|grep -v swpd|awk '{print $13}'").read()
  return int(value)

def GetTotalMemUse():
  value = os.popen("free | grep Mem | awk '{print $3}'").read()
  return int(value) / 1024

def GetVirUse(pid):
  # Mac system will disregard 'h' and print header?
  value = os.popen("ps uh %s | grep %s" % (pid,pid)).read().split()
  if value:
    return int(value[4]) /1024
  return 0

def GetMemUse(pid):
  value = os.popen("ps uh %s | grep %s" % (pid,pid)).read().split()
  if value:
    return int(value[5]) /1024
  return 0

def GetCpuRate(pid):
  value = os.popen("ps uh %s | grep %s" % (pid,pid)).read().split()
  if value:
    return float(value[2])
  return 0

def GetUserName():
  return os.popen("whoami").read().strip()

def GetCondorJobStat(schedd):
  user = GetUserName()
  #return os.popen("condor_q %s -name %s -format \"%%d.\" ClusterId -format \"%%d  \" ProcId -format \"%%d\n\" JobStatus" % (user, shedd)).read().strip()

  cmd = ["condor_q", user, "-name", schedd, "-format", '%d.', "ClusterId", "-format", "%d ", "ProcId", "-format", "%d\n", "JobStatus"]

  #print(cmd)

  output = ""
  for _ in range(100):
    try:
      output = subprocess.check_output(cmd).strip()
      break
    except subprocess.CalledProcessError:
      print("Call '%s' failed. Retry. Waiting 10s."%cmd)
      time.sleep(10)

  #print(output)

  return output

def SubmitCondorJob(script):
  return os.popen("condor_submit %s -verbose | grep \"** Proc\" | awk \'{print $3}\'" % script).read().strip().strip(':')

def GetLSFJobStat():
  user = GetUserName()
  return os.popen("bjobs -u %s" % user).read()

def SubmitLSFJob(script, log, queue, memoryLimit):
  # Job <jobid> is submitted to queue <xxx>.
  #return os.popen("bsub -q %s -o %s %s" % (queue, log, script)).read().split()[1].strip("<>")
  return os.popen("bsub -M %d -q %s %s" % (memoryLimit, queue, script)).read().split()[1].strip("<>")

def GetPBSJobStat():
  user = GetUserName()
  return os.popen("qstat -u %s" % user).read()

def SubmitPBSJob(script, log, walltime):
  # Job <jobid> is submitted to queue <xxx>.                                                                                                                                                                
  #return os.popen("bsub -q %s -o %s %s" % (queue, log, script)).read().split()[1].strip("<>")                                                                                                              
  return os.popen("qsub -l walltime=%s %s" % (walltime, script)).read().split(".")[0].strip("<>")
  #return os.popen("qsub -joe %s" % (script)).read().split(".")[0].strip("<>")

def MakeAndCD(dir):
  if not os.path.isdir(dir):
    try:
      os.makedirs(dir)
    except OSError:
      print('Failed to make dir: ' + dir)
      sys.exit(-1)
  os.chdir(dir)
