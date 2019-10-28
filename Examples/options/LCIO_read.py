#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import LCIODataSvc, CEPCDataSvc

svcname = "LCIODataSvc"
rsvc = LCIODataSvc(svcname, inputs = [
"/cefs/data/FullSim/CEPC240/CEPC_v4/higgs/smart_final_states/E240.Pffh_invi.e0.p0.whizard195//ffh_inv.e0.p0.00001_1000_sim.slcio"
])

wsvc = CEPCDataSvc("EventDataSvc")

from Configurables import PlcioReadAlg
alg = PlcioReadAlg("PlcioReadAlg")
alg.InputCol.Path = "MCParticle"
alg.HeaderCol.Path = "EventHeader"

from Configurables import LCIOInput
lcioinput = LCIOInput("LCIOReader", collections=[
    "EventHeader",
    "MCParticle",
    "TPCCollection"
    #"EventHeader",
    #"MCParticle",
    #"COILCollection",
    #"EcalBarrelSiliconCollection",
    #"EcalBarrelSiliconPreShowerCollection",
    #"EcalEndcapRingCollection",
    #"EcalEndcapRingPreShowerCollection",
    #"EcalEndcapSiliconCollection",
    #"EcalEndcapSiliconPreShowerCollection",
    #"FTD_PIXELCollection",
    #"FTD_STRIPCollection",
    #"HcalBarrelCollection",
    #"HcalEndCapRingsCollection",
    #"HcalEndCapsCollection",
    #"LumiCalCollection",
    #"MuonBarrelCollection",
    #"MuonEndCapCollection",
    #"SETCollection",
    #"SITCollection",
    #"TPCCollection",
    #"TPCSpacePointCollection",
    #"VXDCollection"
    ])
lcioinput.DataSvc = svcname

from Configurables import PodioOutput
plcioout = PodioOutput("PlcioWriter")
plcioout.filename = "lcio2plcio.root"
plcioout.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [lcioinput, alg, plcioout],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [rsvc, wsvc],
                OutputLevel=DEBUG
)
