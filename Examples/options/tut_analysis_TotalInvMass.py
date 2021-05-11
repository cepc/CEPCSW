#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

# read LCIO files
from Configurables import LCIOInput
lcioinput = LCIOInput("LCIOInput")

import glob

# inputlists = ["/cefs/higgs/yudan/CEPC240/Reco_tpc_1800/qqh/Reco_qqh__00001.slcio"]
inputlists = ["/cefs/higgs/yudan/CEPC240/Reco_tpc_1800/qqh/Reco_qqh__00003.slcio"]
# inputlists = glob.glob("/cefs/higgs/yudan/CEPC240/Reco_tpc_1800/qqh/Reco_qqh__*.slcio")

lcioinput.inputs = inputlists
lcioinput.collections = [
    "MCParticle:MCParticle",
    "CalorimeterHit:ECALBarrel",
    "CalorimeterHit:ECALEndcap",
    "CalorimeterHit:HCALBarrel",
    "CalorimeterHit:HCALEndcap",
    "CalorimeterHit:HCALOther",
    "ReconstructedParticle:AncientPFOs",
    "ReconstructedParticle:ArborLICHPFOs"
]

from Configurables import TotalInvMass
total_inv_mass = TotalInvMass("TotalInvMass")

# # write PODIO file
# from Configurables import PodioOutput
# write = PodioOutput("write")
# write.filename = "test.root"
# write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
        TopAlg = [lcioinput, total_inv_mass],
        EvtSel = 'NONE',
        EvtMax = -1,
        ExtSvc = [dsvc],
        OutputLevel=INFO #DEBUG
)

