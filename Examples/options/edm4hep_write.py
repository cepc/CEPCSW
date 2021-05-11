#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import Edm4hepWriteAlg
alg = Edm4hepWriteAlg("Edm4hepWriteAlg")
alg.HeaderOut.Path = "EventHeader"
alg.MCParticleOut.Path = "MCParticle"

from Configurables import PodioOutput
out = PodioOutput("out")
out.filename = "test.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [alg, out],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc=[dsvc],
                OutputLevel=DEBUG
)
