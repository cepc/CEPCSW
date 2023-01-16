from valprod.utils.shellUtil import *
import subprocess

class PSMonitor():
  def __init__(self, name, interval, pid, plotBackend):
    self._sub_monitors = []
    self._sub_monitors.append(VirtMonitor(name,interval,pid,plotBackend))
    self._sub_monitors.append(ResMonitor(name,interval,pid,plotBackend))
    self._sub_monitors.append(CpuMonitor(name,interval,pid,plotBackend))

  def initialize(self):
    for monitor in self._sub_monitors:
      monitor.initialize()

  def execute(self):
    for monitor in self._sub_monitors:
      monitor.execute()

  def finalize(self):
    for monitor in self._sub_monitors:
      monitor.finalize()

class PrmonMonitor():
  def __init__(self, name, interval, pid):
    self._name = name
    self._pid = pid
    self._valid = True
    self._interval = interval

  def initialize(self):
    # Use prmon to monitor the process
    prmon_output = self._name + '.prmon'
    prmon_cmd = ['prmon', '--interval', str(self._interval), '--pid', str(self._pid), '--filename', prmon_output]
    try:
      self.prmon_subprocess = subprocess.Popen(args = prmon_cmd)
    except:
      self._valid = False
      print("prmon execution failed. Switching off monitor")

  def execute(self):
    pass

  def finalize(self):
    if self._valid:
      prmon_rc = self.prmon_subprocess.wait()
      if 0==prmon_rc:
        # Draw figures using the prmon_plot.py tool
        os.system("prmon_plot.py --input %s --xvar wtime --yvar vmem,pss,rss,swap --yunit GB" % prmon_output)
        os.system("prmon_plot.py --input %s --xvar wtime --yvar vmem,pss,rss,swap --diff --yunit MB" % prmon_output)
        os.system("prmon_plot.py --input %s --xvar wtime --yvar utime,stime --yunit SEC --diff --stacked" % prmon_output)
      else:
        print("prmon process failed.")

class PerfMonitor():
  def __init__(self, name, pid):
    self._name = name
    self._perfout = name+'.perfdata'
    self._pid = pid
    self._valid = True

  def initialize(self):
    perf_cmd = ['perf', 'record', '-F 99', '-g', '-p', str(self._pid), '-o', self._perfout]
    try:
      self.perf_subprocess = subprocess.Popen(args = perf_cmd)
    except:
      self._valid = False
      print("perf execution failed. Switching off perf")

  def execute(self):
    pass

  def finalize(self):
    if self._valid:
      perf_rc = self.perf_subprocess.wait()
      if 0==perf_rc:
        ret = os.system("perf script -i %s | stackcollapse-perf.pl > %s.perf-folded" % (self._perfout, self._name))
        if ret:
          print("perf profiling failed.")
          return
        ret = os.system("flamegraph.pl %s.perf-folded > %s_perf.svg" % (self._name, self._name))
        if ret:
          print("perf profiling failed.")
          return
      else:
        print("perf process failed.")
      os.remove(self._perfout)
      os.remove("%s.perf-folded"%self._name)

class PidMonitor():
  def __init__(self, name, interval, pid, backend):
    self.interval = interval
    self.backend = backend
    self.pid = pid
    self.test_name = name
    self.min = 0.
    self.max = 1e15
    self.results = []
    self.hist = None

  def initialize(self):
    pass

  def execute(self):
    sum_value = 0.0

    for pid in self.pid:
      value = eval(self.fun + "(%s)" % pid)
      if value < self.min:
        self.min = value
      if value > self.max:
        self.max = value
      sum_value = sum_value + value
    self.results.append(sum_value)

  def finalize(self):
    if self.backend == 'matplotlib':
      import matplotlib.pyplot as plt
      import numpy as np
      ntime = len(self.results)
      time_seq = np.linspace(0, ntime * self.interval, ntime)
      plt.clf()
      plt.plot(time_seq, self.results)
      plt.xlabel(self.xtitle)
      plt.ylabel(self.ytitle)
      plt.title(self.title)
      plt.savefig(self.test_name + '_' + self.monitor_name + '.png')

    if self.backend == 'root':
      from ROOT import TCanvas, TH1F, TImage
      ntime = len(self.results)
      canvas = TCanvas()
      hist = TH1F(self.test_name, self.title, ntime + 1, 0., ntime * self.interval)
      hist.GetXaxis().SetTitle(self.xtitle)
      hist.GetYaxis().SetTitle(self.ytitle)
      for i in range(ntime):
        hist.SetBinContent(i, self.results[i])
      hist.SetStats(0)
      hist.Draw()
      img = TImage.Create()
      img.FromPad(canvas)
      img.WriteImage(self.test_name + '_' + self.monitor_name + '.png')

class VirtMonitor(PidMonitor):

  def __init__(self, name, interval, pid, plotBackend):
    self.title = "%s Virtual Memory Usage" %name
    self.xtitle = "Time [s]"
    self.ytitle = "Virtual Memory Usage [MB]"
    self.fun = "GetVirUse"
    self.monitor_name = "VirMem"
    PidMonitor.__init__(self, name, interval, pid, plotBackend)
    

class ResMonitor(PidMonitor):

  def __init__(self, name, interval, pid, plotBackend):
    self.title = "%s Resident Memory Usage" %name
    self.xtitle = "Time [s]"
    self.ytitle = "Resident Memory Usage [MB]"
    self.fun = "GetMemUse"
    self.monitor_name = "ResMem"
    PidMonitor.__init__(self, name, interval, pid, plotBackend)


class CpuMonitor(PidMonitor):

  def __init__(self, name, interval, pid, plotBackend):
    self.title = "%s CPU Utilization" % name
    self.xtitle = "Time [s]"
    self.ytitle = "CPU Utilization [/%]"
    self.fun = "GetCpuRate"
    self.monitor_name = "CPURate"
    PidMonitor.__init__(self, name, interval, pid, plotBackend)

class DiskIOMonitor:
  
  def __init__(self, name, interval, pid, plotBackend):
    self._valid = True
    self._name = name
    self._interval = interval
    self._pid = pid
    self._backend = plotBackend
    self._read_count = []
    self._write_count = []
    self._read_bytes = []
    self._write_bytes = []
      
  def initialize(self):
    try:
      import psutil
    except:
      print("Cannot import psutil. Switching off Disk IO Monitor")
      self._valid = False
    try:
      self._process = psutil.Process(self._pid[0])
    except:
      # process not created or already exit
      self._valid = False

  def execute(self):
    if not self._valid:
      return
    try:
      io_counter = self._process.io_counters()
    except:
      # process already existed, or permission error
      return
    self._read_count.append(float(io_counter.read_count))
    self._write_count.append(float(io_counter.write_count))
    self._read_bytes.append(float(io_counter.read_bytes)/1000)   # kb
    self._write_bytes.append(float(io_counter.write_bytes)/1000) # kb

  def finalize(self):
    if not self._valid:
      return

    # Get per-interval data from accumulated data
    for this_result in [self._read_count, self._write_count, self._read_bytes, self._write_bytes]:
      length = len(this_result)
      if length <= 1:  return # too short
      index = length - 1
      while index > 0:
        this_result[index] = this_result[index] - this_result[index-1]
        index = index - 1
      
    if self._backend == 'matplotlib':
      try:
        import matplotlib.pyplot as plt
        import numpy as np
      except:
        print("failed to import matplotlib.")
        return
      ntime = len( self._read_count)
      time_seq = np.linspace(0, ntime * self._interval, ntime)

      # counters
      plt.clf()
      plt.fill_between(time_seq, self._read_count, label="Read", color="red", alpha=0.5)
      plt.fill_between(time_seq, self._read_count, label="Write", color="blue", alpha=0.5)
      plt.legend()
      plt.xlabel("Time [s]")
      plt.ylabel("Count")
      plt.title("%s IO Operations" %self._name)
      plt.savefig(self._name + '_IOCount.png')

      plt.clf()
      plt.fill_between(time_seq, self._read_bytes, label="Read", color="red", alpha=0.5)
      plt.fill_between(time_seq, self._write_bytes, label="Write", color="blue", alpha=0.5)
      plt.legend()
      plt.xlabel("Time [s]")
      plt.ylabel("KB")
      plt.title("%s IO Throughput" %self._name)
      plt.savefig(self._name + '_IOThroughput.png')

    elif self._backend == 'root':
      try:
        from ROOT import TCanvas, TGraph, TImage, TLegend
        from array import array
      except:
        print("failed to import ROOT")
        return
      ntime = len(self._read_count)
      time_seq = list(range(0, ntime * self._interval, self._interval))

      # Decide which
      ts = array('d')
      rc, wc, rb, wb= array('d'), array('d'), array('d'), array('d')
      for i in time_seq:
        ts.append(i)
      for i in self._read_count:
        rc.append(i)
      for i in self._write_count:
        wc.append(i)
      for i in self._read_bytes:
        rb.append(i)
      for i in self._write_bytes:
        wb.append(i)

      canvas = TCanvas()
      graph1 = TGraph(ntime, ts, rc)
      graph2 = TGraph(ntime, ts, wc)
      graph1.SetLineColor(2)
      graph1.SetLineWidth(3)
      graph2.SetLineColor(4)
      graph2.SetLineWidth(3)
      leg = TLegend(.73,.88,.97,.78)
      leg.SetBorderSize(0)
      leg.SetFillColor(0)
      leg.SetFillStyle(0)
      leg.SetTextFont(42)
      leg.SetTextSize(0.035)
      leg.AddEntry(graph1,"Read","L")
      leg.AddEntry(graph2,"Write","L")
      if max(rc) > max(wc):
        graph1.SetTitle("%s IO Operations" %self._name)
        graph1.GetXaxis().SetTitle("Time [s]")
        graph1.GetYaxis().SetTitle("Count")
        graph1.Draw("ACP")
        graph2.Draw("CP")
      else:
        graph2.SetTitle("%s IO Operations" %self._name)
        graph2.GetXaxis().SetTitle("Time [s]")
        graph2.GetYaxis().SetTitle("Count")
        graph2.Draw("ACP")
        graph1.Draw("CP")
      leg.Draw()
      img = TImage.Create()
      img.FromPad(canvas)
      img.WriteImage(self._name + '_IOCount.png')

      canvas = TCanvas()
      graph1 = TGraph(ntime, ts, rb)
      graph2 = TGraph(ntime, ts, wb)
      graph1.SetLineColor(2)
      graph1.SetLineWidth(3)
      graph2.SetLineColor(4)
      graph2.SetLineWidth(3)
      graph1.SetTitle("%s IO Operations" %self._name)
      graph1.GetXaxis().SetTitle("Time [s]")
      graph1.GetYaxis().SetTitle("KB")
      leg = TLegend(.73,.88,.97,.78)
      leg.SetBorderSize(0)
      leg.SetFillColor(0)
      leg.SetFillStyle(0)
      leg.SetTextFont(42)
      leg.SetTextSize(0.035)
      leg.AddEntry(graph1,"Read","L")
      leg.AddEntry(graph2,"Write","L")
      if max(rb) > max(wb):
        graph1.SetTitle("%s IO Operations" %self._name)
        graph1.GetXaxis().SetTitle("Time [s]")
        graph1.GetYaxis().SetTitle("KB")
        graph1.Draw("ACP")
        graph2.Draw("CP")
      else:
        graph2.SetTitle("%s IO Operations" %self._name)
        graph2.GetXaxis().SetTitle("Time [s]")
        graph2.GetYaxis().SetTitle("KB")
        graph2.Draw("ACP")
        graph1.Draw("CP")
      leg.Draw()
      img = TImage.Create()
      img.FromPad(canvas)
      img.WriteImage(self._name + '_IOOperation.png')
    else:
      print("Invalid plotting backend: %s" % self._backend)
