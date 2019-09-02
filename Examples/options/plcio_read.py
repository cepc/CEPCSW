#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import CEPCDataSvc
dsvc = CEPCDataSvc("EventDataSvc", input="test.root")

from Configurables import PlcioReadAlg
alg = PlcioReadAlg("PlcioReadAlg")
alg.InputCol.Path = "MCParticleCol"

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
    "EventHeaderCol",
    "MCParticleCol"
    ])

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
                OutputLevel=DEBUG
)
