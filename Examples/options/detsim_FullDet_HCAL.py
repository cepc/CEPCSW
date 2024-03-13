#!/usr/bin/env python


#######  WARNING #########################################
# 
# This script can not run in CEPCSW master due to geometry
# But can be an example for having truth info in track with FullLDCTrackingAlg
#
###########################################################

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_

seed = [135]
# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_("RndmGenSvc.Engine") # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = seed

rndmgensvc = RndmGenSvc("RndmGenSvc")
rndmgensvc.Engine = rndmengine.name()

#geometry_option = "CRD_o1_v01/CRD_o1_v01.xml"
geometry_option = "CRD_o1_v01/CRD_o1_v01_HCAL.xml"
#...

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

from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
#ntsvc.Output = ["MyTuples DATAFILE='result.root' OPT='NEW' TYP='ROOT'"]

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
#gun.Particles = ["nu_e"]
gun.PositionXs = [0.,0.]
gun.PositionYs = [0.,0.]
gun.PositionZs = [0.,0.]
gun.EnergyMins = [5.,5.] # GeV
gun.EnergyMaxs = [5.,5.] # GeV
gun.ThetaMins  = [90.,91.]   # deg
gun.ThetaMaxs  = [90.,91.]   # deg
gun.PhiMins    = [0.,0.]   # deg
gun.PhiMaxs    = [0.,0.]   # deg


# stdheprdr = StdHepRdr("StdHepRdr")
# stdheprdr.Input = "/cefs/data/stdhep/CEPC240/higgs/exclusive/E240.Pnnh_bb.e0.p0.whizard195/nnh_bb.e0.p0.00001.stdhep"

# lciordr = SLCIORdr("SLCIORdr")
# lciordr.Input = "/cefs/data/stdhep/lcio250/signal/Higgs/E250.Pbbh.whizard195/E250.Pbbh_X.e0.p0.whizard195/Pbbh_X.e0.p0.00001.slcio"
# hepmcrdr = HepMCRdr("HepMCRdr")
# hepmcrdr.Input = "example_UsingIterators.txt"

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtGunTool"]
# genalg.GenTools = ["StdHepRdr"]
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
    # example_anatool.name()
#    "ExampleAnaElemTool",
    "Edm4hepWriterAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"
detsimalg.PhysicsList = "QGSP_BERT_EMV"

from Configurables import MarlinEvtSeeder
evtseeder = MarlinEvtSeeder("EventSeeder")

from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")
#gearsvc.GearXMLFile = "../../Detector/DetCEPCv4/compact/FullDetGear.xml"

from Configurables import TrackSystemSvc
tracksystemsvc = TrackSystemSvc("TrackSystemSvc")

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")

from Configurables import TimeProjectionChamberSensDetTool
tpc_sensdettool = TimeProjectionChamberSensDetTool("TimeProjectionChamberSensDetTool")
tpc_sensdettool.TypeOption = 1


from Configurables import CalorimeterSensDetTool
from Configurables import DriftChamberSensDetTool
cal_sensdettool = CalorimeterSensDetTool("CalorimeterSensDetTool")
cal_sensdettool.CalNamesMergeDisable = ["CaloDetector"]
cal_sensdettool.CalNamesApplyBirks = ["HcalBarrel"]

# digitization
vxdhitname  = "VXDTrackerHits"
sithitname  = "SITTrackerHits"
tpchitname  = "TPCTrackerHits"
sethitname  = "SETTrackerHits"
setspname   = "SETSpacePoints"
ftdhitname  = "FTDTrackerHits"
ftdspname   = "FTDSpacePoints"
from Configurables import PlanarDigiAlg
digiVXD = PlanarDigiAlg("VXDDigi")
digiVXD.SimTrackHitCollection = "VXDCollection"
digiVXD.TrackerHitCollection = vxdhitname
digiVXD.TrackerHitAssociationCollection = "VXDTrackerHitAssociation"
digiVXD.ResolutionU = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.ResolutionV = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.UsePlanarTag = True
#digiVXD.OutputLevel = DEBUG

digiSIT = PlanarDigiAlg("SITDigi")
digiSIT.IsStrip = False
digiSIT.SimTrackHitCollection = "SITCollection"
digiSIT.TrackerHitCollection = sithitname
digiSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
digiSIT.ResolutionU = [0.0072]
digiSIT.ResolutionV = [0.086]
digiSIT.UsePlanarTag = True
#digiSIT.OutputLevel = DEBUG

digiSET = PlanarDigiAlg("SETDigi")
digiSET.IsStrip = False
digiSET.SimTrackHitCollection = "SETCollection"
digiSET.TrackerHitCollection = sethitname
digiSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
digiSET.ResolutionU = [0.0072]
digiSET.ResolutionV = [0.086]
digiSET.UsePlanarTag = True
#digiSET.OutputLevel = DEBUG

# two strip tracker hits -> one space point
from Configurables import SpacePointBuilderAlg
spSET = SpacePointBuilderAlg("SETBuilder")
spSET.TrackerHitCollection = sethitname
spSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
spSET.SpacePointCollection = setspname
spSET.SpacePointAssociationCollection = "SETSpacePointAssociation"
#spSET.OutputLevel = DEBUG


digiFTD = PlanarDigiAlg("FTDDigi")
digiFTD.IsStrip = False
digiFTD.SimTrackHitCollection = "FTDCollection"
digiFTD.TrackerHitCollection = ftdhitname
digiFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
digiFTD.ResolutionU = [0.003, 0.003, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072]
digiFTD.ResolutionV = [0.003, 0.003, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072]
digiFTD.UsePlanarTag = True
#digiFTD.OutputLevel = DEBUG

# two strip tracker hits -> one space point
from Configurables import SpacePointBuilderAlg
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
digiTPC.TPCTrackerHitAssCol = "TPCTrackerHitAssociation"
#digiTPC.OutputLevel = DEBUG

# tracking
from Configurables import SiliconTrackingAlg
tracking = SiliconTrackingAlg("SiliconTracking")
tracking.HeaderCol = "EventHeader"
tracking.VTXHitCollection = vxdhitname
tracking.SITHitCollection = sithitname
tracking.FTDPixelHitCollection = ftdhitname
tracking.FTDSpacePointCollection = ftdspname
tracking.SITRawHitCollection = sithitname
tracking.FTDRawHitCollection = ftdhitname
tracking.UseSIT = True
tracking.SmoothOn = False
tracking.DumpTime = False
tracking.NDivisionsInTheta = 10
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
forward.DumpTime = False
#forward.OutputLevel = DEBUG

from Configurables import TrackSubsetAlg
subset = TrackSubsetAlg("TrackSubset")
subset.TrackInputCollections = ["ForwardTracks", "SiTracks"]
subset.RawTrackerHitCollections = [vxdhitname, sithitname, ftdhitname, ftdspname]
subset.TrackSubsetCollection = "SubsetTracks"
subset.DumpTime = False
#subset.OutputLevel = DEBUG

#TODO: DC reconstruction, as preliminary, use Clupatra like as TPC
from Configurables import ClupatraAlg
clupatra = ClupatraAlg("Clupatra")
clupatra.TPCHitCollection = tpchitname
#clupatra.DistanceCut = 100.
#clupatra.MaxDeltaChi2 = 100.
#clupatra.Chi2Cut = 150.
#clupatra.OutputLevel = DEBUG

from Configurables import FullLDCTrackingAlg
full = FullLDCTrackingAlg("FullTracking")
full.VTXTrackerHits = vxdhitname
full.SITTrackerHits = sithitname
full.TPCTrackerHits = tpchitname  # add TPC or DC tracker hit here, if TPC or DC track is set by full.TPCTracks
full.SETTrackerHits = sethitname
full.FTDPixelTrackerHits = ftdhitname
full.FTDSpacePoints = ftdspname
full.SITRawHits     = sithitname
full.SETRawHits     = sethitname
full.FTDRawHits     = ftdhitname
full.VTXHitRelCol   = "VXDTrackerHitAssociation"
full.SITHitRelCol   = "SITTrackerHitAssociation"
full.SETHitRelCol   = "SETSpacePointAssociation"
full.FTDHitRelCol   = "FTDTrackerHitAssociation"
full.TPCHitRelCol   = "TPCTrackerHitAssociation"
full.TPCTracks = "ClupatraTracks" # add standalone TPC or DC track here
full.SiTracks  = "SubsetTracks"
full.OutputTracks  = "MarlinTrkTracks"
full.DumpTime = False
full.SITHitToTrackDistance = 3.
full.SETHitToTrackDistance = 5.
full.MinChi2ProbForSiliconTracks = 0
#full.OutputLevel = DEBUG

dedxoption = "BetheBlochEquationDedxSimTool"
from Configurables import DriftChamberSensDetTool
dc_sensdettool = DriftChamberSensDetTool("DriftChamberSensDetTool")
dc_sensdettool.DedxSimTool = dedxoption

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

# output
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "CRD_HCal_GamGam_1deg.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSET, digiTPC, tracking, forward, subset, full, out],
    #TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spSET, digiTPC, tracking, forward, subset, out],
    EvtSel = 'NONE',
    EvtMax = 100,
    ExtSvc = [rndmengine, rndmgensvc, dsvc, evtseeder, geosvc, gearsvc, tracksystemsvc],
    #ExtSvc = [rndmengine, rndmgensvc, dsvc, geosvc],
    OutputLevel=INFO
)
