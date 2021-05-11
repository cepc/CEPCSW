#!/usr/bin/env python

from Gaudi.Configuration import *

NTupleSvc().Output = ["MyTuples DATAFILE='sim-rec-trackerEcal.root' OPT='NEW' TYP='ROOT'"]

from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_()
rndmengine.SetSingleton = True
rndmengine.Seeds = [1]

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import MarlinEvtSeeder
evtseeder = MarlinEvtSeeder("EventSeeder")

geometry_option = "CepC_v4-onlyTrackerECAL.xml"

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

from Configurables import GenAlgo
from Configurables import GtGunTool
from Configurables import StdHepRdr
from Configurables import SLCIORdr
from Configurables import HepMCRdr
from Configurables import GenPrinter

gun = GtGunTool("GtGunTool")
gun.Particles = ["e-"]
#gun.EnergyMins = [1]
#gun.EnergyMaxs = [50]
#gun.ThetaMins = [50]
#gun.ThetaMaxs = [130]
#gun.PhiMins = [-90]
#gun.PhiMaxs = [90]
gun.EnergyMins = [10]
gun.EnergyMaxs = [10]
gun.ThetaMins = [90]
gun.ThetaMaxs = [90]
gun.PhiMins = [0]
gun.PhiMaxs = [0]

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtGunTool"]

from Configurables import DetSimSvc
detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg
detsimalg = DetSimAlg("DetSimAlg")
# detsimalg.VisMacs = ["vis.mac"]
detsimalg.RunCmds = [
#    "/physics_lists/factory/addOptical"
]
detsimalg.PhysicsList = "FTFP_BERT"
detsimalg.AnaElems = ["Edm4hepWriterAnaElemTool"]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")

from Configurables import TimeProjectionChamberSensDetTool
tpc_sensdettool = TimeProjectionChamberSensDetTool("TimeProjectionChamberSensDetTool")
tpc_sensdettool.TypeOption = 1

from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")
gearsvc.GearXMLFile = "Detector/DetCEPCv4/compact/FullDetGear.xml"

from Configurables import TrackSystemSvc
tracksystemsvc = TrackSystemSvc("TrackSystemSvc")

vxdhitname  = "VXDTrackerHits"
sithitname  = "SITTrackerHits"
sitspname   = "SITSpacePoints"
tpchitname  = "TPCTrackerHits"
sethitname  = "SETTrackerHits"
setspname   = "SETSpacePoints"
ftdspname   = "FTDSpacePoints"
ftdhitname = "FTDTrackerHits"
from Configurables import PlanarDigiAlg
digiVXD = PlanarDigiAlg("VXDDigi")
digiVXD.SimTrackHitCollection = "VXDCollection"
digiVXD.TrackerHitCollection = vxdhitname
digiVXD.ResolutionU = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.ResolutionV = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]

digiSIT = PlanarDigiAlg("SITDigi")
digiSIT.IsStrip = 1
digiSIT.SimTrackHitCollection = "SITCollection"
digiSIT.TrackerHitCollection = sithitname
digiSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
digiSIT.ResolutionU = [0.007]
digiSIT.ResolutionV = [0.000]

digiSET = PlanarDigiAlg("SETDigi")
digiSET.IsStrip = 1
digiSET.SimTrackHitCollection = "SETCollection"
digiSET.TrackerHitCollection = sethitname
digiSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
digiSET.ResolutionU = [0.007]
digiSET.ResolutionV = [0.000]

digiFTD = PlanarDigiAlg("FTDDigi")
digiFTD.SimTrackHitCollection = "FTDCollection"
digiFTD.TrackerHitCollection = ftdhitname
digiFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
digiFTD.ResolutionU = [0.003, 0.003, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007]
digiFTD.ResolutionV = [0.003, 0.003, 0,     0,     0,     0,     0,     0    ]
#digiFTD.OutputLevel = DEBUG

from Configurables import SpacePointBuilderAlg
spSIT = SpacePointBuilderAlg("SITBuilder")
spSIT.TrackerHitCollection = sithitname
spSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
spSIT.SpacePointCollection = sitspname
spSIT.SpacePointAssociationCollection = "SITSpacePointAssociation"
#spSIT.OutputLevel = DEBUG

spFTD = SpacePointBuilderAlg("FTDBuilder")
spFTD.TrackerHitCollection = ftdhitname
spFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
spFTD.SpacePointCollection = ftdspname
spFTD.SpacePointAssociationCollection = "FTDSpacePointAssociation"
#spFTD.OutputLevel = DEBUG

from Configurables import TPCDigiAlg
digiTPC = TPCDigiAlg("TPCDigi")
digiTPC.TPCCollection = "TPCCollection"
digiTPC.TPCLowPtCollection = "TPCLowPtCollection"
digiTPC.TPCTrackerHitsCol = tpchitname
#digiTPC.OutputLevel = DEBUG

from Configurables import ClupatraAlg
clupatra = ClupatraAlg("Clupatra")
clupatra.TPCHitCollection = tpchitname
#clupatra.OutputLevel = DEBUG

from Configurables import SiliconTrackingAlg
tracking = SiliconTrackingAlg("SiliconTracking")
tracking.HeaderCol = "EventHeader"
tracking.VTXHitCollection = vxdhitname
tracking.SITHitCollection = sitspname
tracking.FTDPixelHitCollection = ftdhitname
tracking.FTDSpacePointCollection = ftdspname
tracking.SITRawHitCollection = sithitname
tracking.FTDRawHitCollection = ftdhitname
tracking.UseSIT = 1
tracking.SmoothOn = 0
#tracking.OutputLevel = DEBUG

from Configurables import ForwardTrackingAlg
forward = ForwardTrackingAlg("ForwardTracking")
forward.FTDPixelHitCollection = ftdhitname
forward.FTDSpacePointCollection = ftdspname
forward.FTDRawHitCollection = ftdhitname
forward.Chi2ProbCut = 0.0
forward.HitsPerTrackMin = 3
forward.BestSubsetFinder = "SubsetSimple"
forward.Criteria = ["Crit2_DeltaPhi","Crit2_StraightTrackRatio","Crit3_3DAngle","Crit3_ChangeRZRatio","Crit3_IPCircleDist","Crit4_3DAngleChange","Crit4_DistToExtrapolation",
                    "Crit2_DeltaRho","Crit2_RZRatio","Crit3_PT"]
forward.CriteriaMin = [0,  0.9,  0,  0.995, 0,  0.8, 0,   20,  1.002, 0.1,      0,   0.99, 0,    0.999, 0,   0.99, 0]  
forward.CriteriaMax = [30, 1.02, 10, 1.015, 20, 1.3, 1.0, 150, 1.08,  99999999, 0.8, 1.01, 0.35, 1.001, 1.5, 1.01, 0.05] 
#forward.OutputLevel = DEBUG

from Configurables import TrackSubsetAlg
subset = TrackSubsetAlg("TrackSubset")
subset.TrackInputCollections = ["ForwardTracks", "SiTracks"]
subset.RawTrackerHitCollections = [vxdhitname, sithitname, ftdhitname, sitspname, ftdspname]
subset.TrackSubsetCollection = "SubsetTracks"
#subset.OutputLevel = DEBUG

from Configurables import FullLDCTrackingAlg
full = FullLDCTrackingAlg("FullTracking")
full.VTXTrackerHits = vxdhitname
full.SITTrackerHits = sitspname
full.TPCTrackerHits = tpchitname
full.SETTrackerHits = setspname
full.FTDPixelTrackerHits = ftdhitname
full.FTDSpacePoints = ftdspname
full.SITRawHits     = sithitname
full.SETRawHits     = sethitname
full.FTDRawHits     = ftdhitname
full.TPCTracks = "ClupatraTracks"
full.SiTracks  = "SubsetTracks"
full.OutputTracks  = "MarlinTrkTracks"
#full.OutputLevel = DEBUG
'''
from Configurables import DumpMCParticleAlg
dumpMC = DumpMCParticleAlg("DumpMC")
dumpMC.MCParticleCollection = "MCParticle"

from Configurables import DumpTrackAlg
dumpFu = DumpTrackAlg("DumpFu")
dumpFu.TrackCollection = "MarlinTrkTracks"
#dumpFu.OutputLevel = DEBUG

dumpCl = DumpTrackAlg("DumpCl")
dumpCl.TrackCollection = "ClupatraTracks"
#dumpCl.OutputLevel = DEBUG

dumpSu = DumpTrackAlg("DumpSu")
dumpSu.TrackCollection = "SubsetTracks"
#dumpSu.OutputLevel = DEBUG

dumpSi = DumpTrackAlg("DumpSi")
dumpSi.TrackCollection = "SiTracks"
#dumpSi.OutputLevel = DEBUG

dumpFo = DumpTrackAlg("DumpFo")
dumpFo.TrackCollection = "ForwardTracks"
#dumpFo.OutputLevel = DEBUG
'''
############################################################
from Configurables import SimHitMergeAlg
simHitMerge = SimHitMergeAlg("SimHitMergeAlg")
simHitMerge.InputCollections=["EcalBarrelCollection", "EcalEndcapsCollection"]
simHitMerge.OutputCollections=["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged"]
############################################################

from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
caloDigi.ReadLCIO = False 
#caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.CalibrECAL = [46.538, 93.0769]
caloDigi.ECALCollections = ["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged"]
caloDigi.ECALReadOutNames= ["EcalBarrelCollection", "EcalEndcapsCollection"]
caloDigi.DigiECALCollection = ["ECALBarrel", "ECALEndcap"]
caloDigi.HCALCollections = []
caloDigi.HCALReadOutNames= []
caloDigi.DigiHCALCollection = []
caloDigi.EventReportEvery = 100
##############################################################################
from Configurables import PandoraPFAlg

pandoralg = PandoraPFAlg("PandoraPFAlg")
pandoralg.debug              = False
pandoralg.use_dd4hep_geo     = True
pandoralg.use_dd4hep_decoder = True
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
        "MCRecoCaloAssociation:MCRecoCaloAssociationCollection" 
        ]
pandoralg.WriteClusterCollection               = "PandoraClusters"              
pandoralg.WriteReconstructedParticleCollection = "PandoraPFOs" 
pandoralg.WriteVertexCollection                = "PandoraPFANewStartVertices"               

pandoralg.PandoraSettingsDefault_xml = "Reconstruction/PFA/Pandora/PandoraSettingsDefault.xml"
#### Do not chage the collection name, only add or remove ###############
pandoralg.TrackCollections      =  ["MarlinTrkTracks"]
pandoralg.ECalCaloHitCollections=  ["ECALBarrel"          , "ECALEndcap"           , "ECALOther"]
pandoralg.ECalReadOutNames      =  ["EcalBarrelCollection", "EcalEndcapsCollection", "ECALOther"]
pandoralg.HCalCaloHitCollections=  ["HCALBarrel", "HCALEndcap", "HCALOther"]
pandoralg.HCalReadOutNames      =  ["HcalBarrelCollection", "HcalEndcapsCollection", "HCALOther"]
pandoralg.LCalCaloHitCollections=  ["LCAL"]
pandoralg.LCalReadOutNames      =  ["LcalCollection"]
pandoralg.LHCalCaloHitCollections= ["LHCAL"]
pandoralg.LHCalReadOutNames      = ["LHcalCollection"]
pandoralg.MuonCaloHitCollections=  ["MUON"]
pandoralg.MuonCalReadOutNames    = ["MuoncalCollection"]
pandoralg.MCParticleCollections =  ["MCParticle"]
pandoralg.RelCaloHitCollections =  ["MCRecoCaloAssociationCollection"]
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
pandoralg.AbsorberRadLengthECal= 0.2854 #= 1/3.504 mm 
pandoralg.AbsorberIntLengthECal= 0.0101 #= 1/99.46 mm 
pandoralg.AbsorberRadLengthHCal= 0.0569 
pandoralg.AbsorberIntLengthHCal= 0.006  
pandoralg.AbsorberRadLengthOther= 0.0569
pandoralg.AbsorberIntLengthOther= 0.006 

##############################################################################

# write PODIO file
from Configurables import PodioOutput
write = PodioOutput("write")
write.filename = "test.root"
write.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    #TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSIT, spFTD, digiTPC, clupatra, tracking, forward, subset, full, dumpMC, dumpFu, dumpCl, dumpSu, dumpSi, dumpFo, write],
    TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSIT, spFTD, digiTPC, clupatra, tracking, forward, subset, full, simHitMerge, caloDigi, pandoralg, write],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [rndmengine, dsvc, evtseeder, gearsvc, geosvc, tracksystemsvc],
    HistogramPersistency='ROOT',
    OutputLevel=INFO
)
