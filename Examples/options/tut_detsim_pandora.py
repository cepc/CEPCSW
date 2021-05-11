#!/usr/bin/env python

import os
print(os.environ["DD4HEP_LIBRARY_PATH"])
import sys

from Gaudi.Configuration import *

##############################################################################
# Random Number Svc
##############################################################################
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_

seed = [42]

# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_("RndmGenSvc.Engine") # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = seed

rndmgensvc = RndmGenSvc("RndmGenSvc")
rndmgensvc.Engine = rndmengine.name()

##############################################################################
# Event Data Svc
##############################################################################
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")


##############################################################################
# Geometry Svc
##############################################################################

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

##############################################################################
# Physics Generator
##############################################################################
from Configurables import GenAlgo
from Configurables import GtGunTool
from Configurables import StdHepRdr
from Configurables import SLCIORdr
from Configurables import HepMCRdr
from Configurables import GenPrinter

gun = GtGunTool("GtGunTool")
gun.Particles = ["gamma"]
gun.EnergyMins= [10] # GeV
gun.EnergyMaxs= [10] # GeV
gun.ThetaMins = [80] # degree
gun.ThetaMaxs = [80] # degree
gun.PhiMins   = [0 ] # degree
gun.PhiMaxs   = [0 ] # degree


#stdheprdr = StdHepRdr("StdHepRdr")
#stdheprdr.Input = "/cefs/data/stdhep/CEPC250/higgs/E250.Pbbh.whizard195/E250.Pbbh_X.e0.p0.whizard195/Pbbh_X.e0.p0.00001.stdhep"

# lciordr = SLCIORdr("SLCIORdr")
# lciordr.Input = "/cefs/data/stdhep/lcio250/signal/Higgs/E250.Pbbh.whizard195/E250.Pbbh_X.e0.p0.whizard195/Pbbh_X.e0.p0.00001.slcio"

# hepmcrdr = HepMCRdr("HepMCRdr")
# hepmcrdr.Input = "example_UsingIterators.txt"

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtGunTool"]
#genalg.GenTools = ["StdHepRdr"]
# genalg.GenTools = ["StdHepRdr", "GenPrinter"]
# genalg.GenTools = ["SLCIORdr", "GenPrinter"]
# genalg.GenTools = ["HepMCRdr", "GenPrinter"]

##############################################################################
# Detector Simulation
##############################################################################
from Configurables import DetSimSvc

detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg

detsimalg = DetSimAlg("DetSimAlg")
detsimalg.RandomSeeds = seed

# detsimalg.VisMacs = ["vis.mac"]

detsimalg.RunCmds = [
#    "/tracking/verbose 1",
]
detsimalg.AnaElems = [
    "Edm4hepWriterAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")

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
from Configurables import GearSvc
gearSvc  = GearSvc("GearSvc")
gearSvc.GearXMLFile = "Detector/DetCEPCv4/compact/FullDetGear.xml"
##############################################################################
from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
ntsvc.Output = ["MyTuples DATAFILE='detsim_Pan_ana.root' OPT='NEW' TYP='ROOT'"]
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
pandoralg.ECalCaloHitCollections=  ["ECALBarrel", "ECALEndcap", "ECALOther"]
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
        TopAlg = [genalg, detsimalg, simHitMerge, caloDigi, pandoralg, write],
        EvtSel = 'NONE',
        EvtMax = 10,
        ExtSvc = [rndmengine, rndmgensvc, dsvc, geosvc, gearSvc,detsimsvc],
        HistogramPersistency = "ROOT",
        OutputLevel=INFO
)
