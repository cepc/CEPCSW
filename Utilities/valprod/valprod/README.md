## Use valprod to build unit test cases

### Via the jMonitor command

The jMonitor command monitors the executable by adding a thin shell on top of it.
Example:

`jMonitor runSimulation.exe`

The detailed usage:
```
usage: jMonitor [-h] [-T MAXTIME] [-WT WALLTIME] [-M MAXVIR] [--enable-parser]
                [-p FATALPATTERN] [--enable-monitor] [-i INTERVAL]
                [-b {root,matplotlib}] [--enable-plotref] [-f PLOTREFFILE]
                [-o PLOTREFOUTPUT] [-m {Kolmogorov,Chi2}] [-c HISTTESTCUT]
                [--gen-log] [-n NAME]
                command

positional arguments:
  command               job to be monitored

optional arguments:
  -h, --help            show this help message and exit
  -T MAXTIME, --max-time MAXTIME
                        Time limit of process (s), if exceeded, process will
                        be killed
  -WT WALLTIME, --walltime WALLTIME
                        Time limit of process (hh:mm:ss) for PBS, if exceeded,
                        process will be killed
  -M MAXVIR, --max-memory MAXVIR
                        Memory limit of process (Mb), if exceeded, process
                        will be killed
  --enable-parser       If enabled, log of the process will be parsed
  -p FATALPATTERN, --pattern FATALPATTERN
                        Python re patterns for the parser
  --enable-monitor      If enabled, process will be monitored (memory, cpu)
  -i INTERVAL, --interval INTERVAL
                        Time interval for monitor
  -b {root,matplotlib}, --monitor-backend {root,matplotlib}
                        Backend for drawing monitoring figures
  --enable-plotref      If enabled, results of the process will be compared
  -f PLOTREFFILE, --plotref-files PLOTREFFILE
                        reference file for plot testing
  -o PLOTREFOUTPUT, --plotref-output PLOTREFOUTPUT
                        output root file for plot comparison
  -m {Kolmogorov,Chi2}, --histtest-method {Kolmogorov,Chi2}
                        Method of histogram testing
  -c HISTTESTCUT, --histtest-cut HISTTESTCUT
                        P-Value cut for histogram testing
  --gen-log             whether to generate log file
  -n NAME, --name NAME  name of the job
```

### Via the API 

The UnitTest module.   

Example: 
```
       mytest = UnitTest()
       mytest.addCase('test1', 'python run.py')
       mytest.run()
```

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
    options of one single test case, example: `mytest.addCase('test2', 'python run.py', timeLimit=900, genLog=True, RESMonitor=True)`

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
