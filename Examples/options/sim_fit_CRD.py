#!/usr/bin/env python
from Gaudi.Configuration import *

##############################################################################
# Event service
##############################################################################
from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

##############################################################################
# Random seed
##############################################################################
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_() # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = [10]

from Configurables import MarlinEvtSeeder
evtseeder = MarlinEvtSeeder("EventSeeder")

##############################################################################
# Detector geometry
##############################################################################
geometry_option = "CRD_o1_v01/CRD_o1_v01_noECAL.xml"

if not os.getenv("DETCRDROOT"):
    print("Can't find the geometry. Please setup envvar DETCRDROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCRDROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path

from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")

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
gun.Particles = ["mu-"]
gun.EnergyMins = [100.] # GeV
gun.EnergyMaxs = [100.] # GeV
gun.ThetaMins  = [0]    # deg
gun.ThetaMaxs  = [180]  # deg
gun.PhiMins    = [0]    # deg
gun.PhiMaxs    = [360]  # deg
# stdheprdr = StdHepRdr("StdHepRdr")
# stdheprdr.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"
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
#detsimalg.VisMacs = ["vis.mac"]
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

from Configurables import CalorimeterSensDetTool
from Configurables import DriftChamberSensDetTool

calo_sensdettool = CalorimeterSensDetTool("CalorimeterSensDetTool")
driftchamber_sensdettool = DriftChamberSensDetTool("DriftChamberSensDetTool")

##############################################################################
# DedxSimTool
##############################################################################
# dedxoption = "DummyDedxSimTool"
dedxoption = "BetheBlochEquationDedxSimTool"

driftchamber_sensdettool.DedxSimTool = dedxoption

from Configurables import DummyDedxSimTool
from Configurables import BetheBlochEquationDedxSimTool

if dedxoption == "DummyDedxSimTool":
    dedx_simtool = DummyDedxSimTool("DummyDedxSimTool")
elif dedxoption == "BetheBlochEquationDedxSimTool":
    dedx_simtool = BetheBlochEquationDedxSimTool("BetheBlochEquationDedxSimTool")
    dedx_simtool.material_Z = 2
    dedx_simtool.material_A = 4
    dedx_simtool.scale = 10
    dedx_simtool.resolution = 0.0001

##############################################################################
# DCHDigiAlg
##############################################################################
from Configurables import DCHDigiAlg
dCHDigiAlg = DCHDigiAlg("DCHDigiAlg")
dCHDigiAlg.readout = "DriftChamberHitsCollection"
dCHDigiAlg.drift_velocity = 40#um/ns
dCHDigiAlg.mom_threshold = 0 #GeV
dCHDigiAlg.WriteAna  = True

##############################################################################
# PODIO
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "CRD-o1-v01-sim00.root"
out.outputCommands = ["keep *"]


##############################################################################
# Silicon digitization
##############################################################################
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


##############################################################################
# Silicon reconstruction
##############################################################################
from Configurables import TrackSystemSvc
tracksystemsvc = TrackSystemSvc("TrackSystemSvc")

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

##############################################################################
# TruthTrackerAlg
##############################################################################
from Configurables import TruthTrackerAlg
truthTrackerAlg = TruthTrackerAlg("TruthTrackerAlg")
truthTrackerAlg.SiSubsetTrackCollection = "SubsetTracks"

##############################################################################
# RecGenfitAlgSDT
##############################################################################
from Configurables import RecGenfitAlgSDT
recGenfitAlgSDT = RecGenfitAlgSDT("RecGenfitAlgSDT")
recGenfitAlgSDT.debug=10

##############################################################################
# NTupleSvc
##############################################################################
from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
ntsvc.Output = [
#"MyTuples DATAFILE='DCH_digi_ana.root' OPT='NEW' TYP='ROOT'",
                "RecGenfitAlgSDT DATAFILE='fit_SDT.root' OPT='NEW' TYP='ROOT'"]


##############################################################################
# ApplicationMgr
##############################################################################
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSIT,
    spFTD, tracking, forward, subset, dCHDigiAlg, truthTrackerAlg,
    recGenfitAlgSDT, out],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [rndmengine, dsvc, ntsvc, evtseeder, geosvc, gearsvc, tracksystemsvc],
    HistogramPersistency = "ROOT",
    OutputLevel=DEBUG
)
