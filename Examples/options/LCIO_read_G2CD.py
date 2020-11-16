#!/usr/bin/env python

from Gaudi.Configuration import *
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

# read LCIO files
from Configurables import LCIOInput
read = LCIOInput("read")
read.inputs = [
"/cefs/higgs/wxfang/cepc/FS/FS_Barrel_Mom150_1000MeV_x_theta_phi_mom_bin_resVsMom_Etuning/sim_111_0830.slcio"
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
#########################################################################
geometry_option = "CepC_v4-onlyECAL.xml"

if not os.getenv("DETCEPCV4ROOT"):
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path
########################################################################
from Configurables import GearSvc
gearSvc  = GearSvc("GearSvc")
gearSvc.GearXMLFile = "Detector/DetCEPCv4/compact/FullDetGear.xml"
##############################################################################
from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
caloDigi.ReadLCIO = True 
caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.ECALCollections = ["EcalBarrelSiliconCollection"]
caloDigi.DigiECALCollection = ["EcalBarrel"]
caloDigi.HCALCollections = []
caloDigi.EventReportEvery = 1

##############################################################################

# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "test.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
        #TopAlg = [read, caloDigi, write],
        TopAlg = [read, caloDigi],
        EvtSel = 'NONE',
        EvtMax = 10,
        ExtSvc = [dsvc, geosvc, gearSvc],
        OutputLevel=INFO
)
