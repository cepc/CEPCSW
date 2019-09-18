#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import LCIODataSvc
dsvc = LCIODataSvc("EventDataSvc", input="/cefs/data/FullSim/CEPC240/CEPC_v4/higgs/E240.Pe2e2h_bb.e0.p0.whizard195/e2e2h_bb.e0.p0.00001_000000_sim.slcio")

from Configurables import PlcioReadAlg
alg = PlcioReadAlg("PlcioReadAlg")
alg.InputCol.Path = "MCParticle"
alg.HeaderCol.Path = "EventHeader"

from Configurables import LCIOInput
lcioinput = LCIOInput("LCIOReader", collections=[
    "EventHeader",
    "MCParticle",
    "COILCollection",
    "EcalBarrelSiliconCollection",
    "EcalBarrelSiliconPreShowerCollection",
    "EcalEndcapRingCollection",
    "EcalEndcapRingPreShowerCollection",
    "EcalEndcapSiliconCollection",
    "EcalEndcapSiliconPreShowerCollection",
    "FTD_PIXELCollection",
    "FTD_STRIPCollection",
    "HcalBarrelCollection",
    "HcalEndCapRingsCollection",
    "HcalEndCapsCollection",
    "LumiCalCollection",
    "MuonBarrelCollection",
    "MuonEndCapCollection",
    "SETCollection",
    "SITCollection",
    "TPCCollection",
    "TPCSpacePointCollection",
    "VXDCollection"
    ])

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [lcioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
                OutputLevel=DEBUG
)
