#!/usr/bin/env python

from Gaudi.Configuration import *
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

# read LCIO files
from Configurables import LCIOInput
read = LCIOInput("read")
read.inputs = [
#"/cefs/higgs/wxfang/cepc/Pandora/ele_lcio_full_det_rec/reco_sim_0.slcio"
"/cefs/higgs/wxfang/cepc/Pandora/ele_lcio_full_det_rec/reco_sim_34.slcio"
]
read.collections = [
        "MCParticle:MCParticle",
        #"SimCalorimeterHit:EcalBarrelSiliconCollection",
        "CalorimeterHit:ECALBarrel",
        "CalorimeterHit:ECALEndcap",
        "CalorimeterHit:ECALOther" ,
        ######### HCAL will effect the reco efficiency close to gap region ######
        "CalorimeterHit:HCALBarrel",
        "CalorimeterHit:HCALEndcap",
        "CalorimeterHit:HCALOther",
        #"TrackerHit:VXDTrackerHits",
        #"TrackerHit:SITTrackerHits",
        "TrackerHit:SITSpacePoints",
        "TrackerHit:TPCTrackerHits",
        #"TrackerHit:SETTrackerHits",
        "TrackerHit:SETSpacePoints",
        #"TrackerHit:FTDStripTrackerHits",
        "TrackerHit:FTDSpacePoints",
        #"TrackerHit:FTDPixelTrackerHits",
        "Track:ClupatraTrackSegments", 
        "Track:ClupatraTracks", 
        "Track:ForwardTracks", 
        "Track:SiTracks", 
        "Track:SubsetTracks",
        "Track:MarlinTrkTracks", 
        "Vertex:KinkVertices",
        "Vertex:ProngVertices",
        "Vertex:V0Vertices",
        "ReconstructedParticle:KinkRecoParticles",
        "ReconstructedParticle:ProngRecoParticles",
        "ReconstructedParticle:V0RecoParticles"
]
##############################################################################
from Configurables import GearSvc
gearSvc  = GearSvc("GearSvc")
gearSvc.GearXMLFile = "Detector/DetCEPCv4/compact/FullDetGear.xml"
##############################################################################
from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
ntsvc.Output = ["MyTuples DATAFILE='LCIO_Pan_ana.root' OPT='NEW' TYP='ROOT'"]
##############################################################################
from Configurables import PandoraPFAlg

pandoralg = PandoraPFAlg("PandoraPFAlg")
pandoralg.debug              = False
pandoralg.use_dd4hep_geo     = False
pandoralg.use_dd4hep_decoder = False
pandoralg.use_preshower      = False
pandoralg.WriteAna           = True
pandoralg.collections = [
        "MCParticle:MCParticle",
        "CalorimeterHit:ECALBarrel",
        "CalorimeterHit:ECALEndcap",
        "CalorimeterHit:ECALOther" ,
        "CalorimeterHit:HCALBarrel",
        "CalorimeterHit:HCALEndcap",
        "CalorimeterHit:HCALOther" ,
        "CalorimeterHit:MUON", 
        "CalorimeterHit:LCAL", 
        "CalorimeterHit:LHCAL", 
        "CalorimeterHit:BCAL", 
        "Vertex:KinkVertices", 
        "Vertex:ProngVertices", 
        "Vertex:SplitVertices", 
        "Vertex:V0Vertices", 
        "Track:MarlinTrkTracks", 
        "MCRecoCaloAssociation:RecoCaloAssociation_ECALBarrel" 
        ]
pandoralg.WriteClusterCollection               = "PandoraClusters"              
pandoralg.WriteReconstructedParticleCollection = "PandoraPFOs" 
pandoralg.WriteVertexCollection                = "PandoraPFANewStartVertices"               

pandoralg.PandoraSettingsDefault_xml = "Reconstruction/PFA/Pandora/PandoraSettingsDefault.xml"
#### Do not chage the collection name, only add or delete ###############
pandoralg.TrackCollections      =  ["MarlinTrkTracks"]
pandoralg.ECalCaloHitCollections=  ["ECALBarrel", "ECALEndcap", "ECALOther"]
pandoralg.HCalCaloHitCollections=  ["HCALBarrel", "HCALEndcap", "HCALOther"]
pandoralg.LCalCaloHitCollections=  ["LCAL"]
pandoralg.LHCalCaloHitCollections= ["LHCAL"]
pandoralg.MuonCaloHitCollections=  ["MUON"]
pandoralg.MCParticleCollections =  ["MCParticle"]
pandoralg.RelCaloHitCollections =  ["RecoCaloAssociation_ECALBarrel", "RecoCaloAssociation_ECALEndcap", "RecoCaloAssociation_ECALOther", "RecoCaloAssociation_HCALBarrel", "RecoCaloAssociation_HCALEndcap", "RecoCaloAssociation_HCALOther", "RecoCaloAssociation_LCAL", "RecoCaloAssociation_LHCAL", "RecoCaloAssociation_MUON"]
pandoralg.RelTrackCollections   =  ["MarlinTrkTracksMCTruthLink"]
pandoralg.KinkVertexCollections =  ["KinkVertices"]
pandoralg.ProngVertexCollections=  ["ProngVertices"]
pandoralg.SplitVertexCollections=  ["SplitVertices"]
pandoralg.V0VertexCollections   =  ["V0Vertices"]
pandoralg.ECalToMipCalibration  = 160.0 
pandoralg.HCalToMipCalibration  = 34.8 
pandoralg.ECalMipThreshold      = 0.5 
pandoralg.HCalMipThreshold      = 0.3 
pandoralg.ECalToEMGeVCalibration= 0.9 #for G2CD Digi, 1.007 for NewLDCaloDigi 
pandoralg.HCalToEMGeVCalibration= 1.007 
pandoralg.ECalToHadGeVCalibrationBarrel= 1.12 #very small effect 
pandoralg.ECalToHadGeVCalibrationEndCap= 1.12 
pandoralg.HCalToHadGeVCalibration= 1.07
pandoralg.MuonToMipCalibration= 10.0 
pandoralg.DigitalMuonHits= 0 
pandoralg.MaxHCalHitHadronicEnergy   = 1.0 
pandoralg.UseOldTrackStateCalculation= 0 
pandoralg.AbsorberRadLengthECal= 0.2854 
pandoralg.AbsorberIntLengthECal= 0.0101 
pandoralg.AbsorberRadLengthHCal= 0.0569 
pandoralg.AbsorberIntLengthHCal= 0.006  
pandoralg.AbsorberRadLengthOther= 0.0569
pandoralg.AbsorberIntLengthOther= 0.006 

##############################################################################

# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "pan_test.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
        TopAlg = [read, pandoralg],
        EvtSel = 'NONE',
        EvtMax = 10,
        ExtSvc = [dsvc, gearSvc],
        HistogramPersistency = "ROOT",
        OutputLevel=INFO
)
