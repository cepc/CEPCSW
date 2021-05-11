#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc",
                 # input="test.root"
                 input="test-detsim10.root"
)

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
    "EventHeader",
    "MCParticle",
    "SimCalorimeterHit"
    ])

from Configurables import Edm4hepReadAlg
alg = Edm4hepReadAlg("Edm4hepReadAlg")
#alg.HeaderCol.Path = "EventHeader"
#alg.MCParticleCol.Path = "MCParticle"
alg.SimCalorimeterHitCol.Path = "SimCalorimeterHit"

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
                OutputLevel=DEBUG
)
