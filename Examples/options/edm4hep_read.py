#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import K4DataSvc
dsvc = K4DataSvc("EventDataSvc", input="test.root")

from Configurables import Edm4hepReadAlg
alg = Edm4hepReadAlg("Edm4hepReadAlg")
alg.HeaderCol.Path = "EventHeader"
alg.InputCol.Path = "MCParticle"

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
    "EventHeader",
    "MCParticle"
    ])

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
                OutputLevel=DEBUG
)
