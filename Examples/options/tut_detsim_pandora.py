#!/usr/bin/env python

import os
print(os.environ["DD4HEP_LIBRARY_PATH"])
import sys
# sys.exit(0)

from Gaudi.Configuration import *

##############################################################################
# Random Number Svc
##############################################################################
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_

# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_() # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = [42]

# rndmgensvc = RndmGenSvc("RndmGenSvc")
# rndmgensvc.Engine = rndmengine.name()


##############################################################################
# Event Data Svc
##############################################################################
#from Configurables import CEPCDataSvc
#dsvc = CEPCDataSvc("EventDataSvc")
from Configurables import K4DataSvc
dsvc = K4DataSvc("EventDataSvc")


##############################################################################
# Geometry Svc
##############################################################################

# geometry_option = "CepC_v4-onlyTracker.xml"
geometry_option = "CepC_v4-onlyVXD.xml"

if not os.getenv("DETCEPCV4ROOT"):
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeoSvc
geosvc = GeoSvc("GeoSvc")
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
gun.Particles = ["gamma","gamma"]
gun.EnergyMins= [10, 10] # GeV
gun.EnergyMaxs= [10, 10] # GeV
gun.ThetaMins = [90, 90] # degree
gun.ThetaMaxs = [90, 90] # degree
gun.PhiMins   = [0,  1 ] # degree
gun.PhiMaxs   = [0,  1 ] # degree


stdheprdr = StdHepRdr("StdHepRdr")
#stdheprdr.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"
#stdheprdr.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"
stdheprdr.Input = "/cefs/data/stdhep/CEPC250/higgs/E250.Pbbh.whizard195/E250.Pbbh_X.e0.p0.whizard195/Pbbh_X.e0.p0.00001.stdhep"

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

# from Configurables import ExampleAnaElemTool
# example_anatool = ExampleAnaElemTool("ExampleAnaElemTool")

from Configurables import DetSimAlg

detsimalg = DetSimAlg("DetSimAlg")

# detsimalg.VisMacs = ["vis.mac"]

detsimalg.RunCmds = [
#    "/tracking/verbose 1",
]
detsimalg.AnaElems = [
    # example_anatool.name()
    # "ExampleAnaElemTool"
    "Edm4hepWriterAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")

##############################################################################
from Configurables import CaloDigiAlg
example_CaloDigiAlg = CaloDigiAlg("CaloDigiAlg")
example_CaloDigiAlg.Scale = 1
example_CaloDigiAlg.SimCaloHitCollection = "SimCalorimeterCol"
example_CaloDigiAlg.CaloHitCollection    = "ECALBarrel"
example_CaloDigiAlg.CaloAssociationCollection    = "RecoCaloAssociation_ECALBarrel"
##############################################################################
from Configurables import GearSvc
gearSvc  = GearSvc("GearSvc")
gearSvc.GearXMLFile = "../Detector/DetCEPCv4/compact/FullDetGear.xml"
##############################################################################
from Configurables import PandoraPFAlg

pandoralg = PandoraPFAlg("PandoraPFAlg")
## KEEP same with lcioinput name for the ReadXXX ###########
pandoralg.ReadMCParticle                       = "MCParticle"                   
pandoralg.ReadECALBarrel                       = "ECALBarrel"                   
pandoralg.ReadECALEndcap                       = "ECALEndcap"                   
pandoralg.ReadECALOther                        = "ECALOther"                    
pandoralg.ReadHCALBarrel                       = "HCALBarrel"                   
pandoralg.ReadHCALEndcap                       = "HCALEndcap"                   
pandoralg.ReadHCALOther                        = "HCALOther"                    
pandoralg.ReadMUON                             = "MUON"                         
pandoralg.ReadLCAL                             = "LCAL"                         
pandoralg.ReadLHCAL                            = "LHCAL"                        
pandoralg.ReadBCAL                             = "BCAL"                         
pandoralg.ReadKinkVertices                     = "KinkVertices"                 
pandoralg.ReadProngVertices                    = "ProngVertices"                
pandoralg.ReadSplitVertices                    = "SplitVertices"                
pandoralg.ReadV0Vertices                       = "V0Vertices"                   
pandoralg.ReadTracks                           = "MarlinTrkTracks"                       
pandoralg.MCRecoCaloAssociation                = "RecoCaloAssociation_ECALBarrel"                       
pandoralg.WriteClusterCollection               = "PandoraClusters"              
pandoralg.WriteReconstructedParticleCollection = "PandoraPFOs" 
pandoralg.WriteVertexCollection                = "PandoraPFANewStartVertices"               
pandoralg.AnaOutput = "Ana.root"

pandoralg.PandoraSettingsDefault_xml = "../Reconstruction/PFA/Pandora/PandoraSettingsDefault.xml"
#### Do not chage the collection name, only add or remove ###############
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
        TopAlg = [genalg, detsimalg, example_CaloDigiAlg, pandoralg],
        EvtSel = 'NONE',
        EvtMax = 10,
        ExtSvc = [rndmengine, dsvc, geosvc, gearSvc,detsimsvc],
        OutputLevel=INFO
)
