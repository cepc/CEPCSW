#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import TBDataClusterReader
alg = TBDataClusterReader("TBDataClusterReader")
alg.TrackerHitOut.Path = "TrackerHit"

from Configurables import PodioOutput
out = PodioOutput("out")
out.filename = "Testbeamcluster.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [alg, out],
                EvtSel = 'NONE',
                EvtMax = 1883792,
                ExtSvc=[dsvc],
                # OutputLevel=DEBUG
)
