#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import TBDataReader
alg = TBDataReader("TBDataReader")
alg.RawTimeSeriesOut.Path = "RawTimeSeries"

from Configurables import PodioOutput
out = PodioOutput("out")
out.filename = "Testbeamdata.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [alg, out],
                EvtSel = 'NONE',
                EvtMax = 8285,
                ExtSvc=[dsvc],
                # OutputLevel=DEBUG
)
