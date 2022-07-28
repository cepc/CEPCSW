#!/usr/bin/env python

from Gaudi.Configuration import *
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

# read LCIO files
from Configurables import LCIOInput
read = LCIOInput("read")
read.inputs = [
"/cefs/data/FullSim/CEPC240/CEPC_v4/higgs/smart_final_states/E240.Pffh_invi.e0.p0.whizard195/ffh_inv.e0.p0.00001_1000_sim.slcio"
]
read.collections = [
        "MCParticle:MCParticle",
        "SimCalorimeterHit:EcalBarrelSiliconCollection",
        #"CalorimeterHit:ECALBarrel",
        #"CalorimeterHit:ECALEndcap",
        #"CalorimeterHit:ECALOther" ,
        ########## HCAL will effect the reco efficiency close to gap region ######
        #"CalorimeterHit:HCALBarrel",
        #"CalorimeterHit:HCALEndcap",
        #"CalorimeterHit:HCALOther",
        ##"TrackerHit:VXDTrackerHits",
        ##"TrackerHit:SITTrackerHits",
        #"TrackerHit:SITSpacePoints",
        #"TrackerHit:TPCTrackerHits",
        ##"TrackerHit:SETTrackerHits",
        #"TrackerHit:SETSpacePoints",
        ##"TrackerHit:FTDStripTrackerHits",
        #"TrackerHit:FTDSpacePoints",
        ##"TrackerHit:FTDPixelTrackerHits",
        #"Track:ClupatraTrackSegments", 
        #"Track:ClupatraTracks", 
        #"Track:ForwardTracks", 
        #"Track:SiTracks", 
        #"Track:SubsetTracks",
        #"Track:MarlinTrkTracks", 
        #"Vertex:KinkVertices",
        #"Vertex:ProngVertices",
        #"Vertex:V0Vertices",
        #"ReconstructedParticle:KinkRecoParticles",
        #"ReconstructedParticle:ProngRecoParticles",
        #"ReconstructedParticle:V0RecoParticles"
]

from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "lcio2plcio.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [read, write],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
                OutputLevel=DEBUG
)
