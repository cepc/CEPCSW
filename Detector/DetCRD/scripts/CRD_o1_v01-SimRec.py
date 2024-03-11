#!/usr/bin/env python
import os

from Gaudi.Configuration import *

NTupleSvc().Output = ["MyTuples DATAFILE='test.root' OPT='NEW' TYP='ROOT'"]

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
seed = [10]
# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_("RndmGenSvc.Engine") # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = seed

rndmgensvc = RndmGenSvc("RndmGenSvc")
rndmgensvc.Engine = rndmengine.name()

# option for standalone tracker study 
geometry_option = "CRD_o1_v01/CRD_o1_v01-onlyTracker.xml"
#geometry_option = "CRD_o1_v01/CRD_o1_v01.xml"

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
gun.ThetaMins  = [85]    # deg
gun.ThetaMaxs  = [85]  # deg
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
detsimalg.RandomSeeds = seed
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
    
from Configurables import MarlinEvtSeeder
evtseeder = MarlinEvtSeeder("EventSeeder")

from Configurables import GearSvc
gearsvc = GearSvc("GearSvc")

from Configurables import TrackSystemSvc
tracksystemsvc = TrackSystemSvc("TrackSystemSvc")

# digitization
vxdhitname  = "VXDTrackerHits"
sithitname  = "SITTrackerHits"
dchitname   = "DCTrackerHits"
sethitname  = "SETTrackerHits"
setspname   = "SETSpacePoints"
ftdhitname  = "FTDTrackerHits"
ftdspname   = "FTDSpacePoints"
from Configurables import PlanarDigiAlg
digiVXD = PlanarDigiAlg("VXDDigi")
digiVXD.SimTrackHitCollection = "VXDCollection"
digiVXD.TrackerHitCollection = vxdhitname
digiVXD.ResolutionU = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.ResolutionV = [0.0028, 0.006, 0.004, 0.004, 0.004, 0.004]
digiVXD.UsePlanarTag = True
# if Parameterized spatial resolution, set ParameterizeResolution to True
digiVXD.ParameterizeResolution = False
digiVXD.ParametersU = [5.60959e-03, 5.74913e-03, 7.03433e-03, 1.99516, -663.952, 3.752e-03, 0, -0.0704734, 0.0454867e-03, 1.07359]
digiVXD.ParametersV = [5.60959e-03, 5.74913e-03, 7.03433e-03, 1.99516, -663.952, 3.752e-03, 0, -0.0704734, 0.0454867e-03, 1.07359]
#digiVXD.OutputLevel = DEBUG

digiSIT = PlanarDigiAlg("SITDigi")
digiSIT.IsStrip = False
digiSIT.SimTrackHitCollection = "SITCollection"
digiSIT.TrackerHitCollection = sithitname
digiSIT.TrackerHitAssociationCollection = "SITTrackerHitAssociation"
digiSIT.ResolutionU = [0.0072]
digiSIT.ResolutionV = [0.086]
digiSIT.UsePlanarTag = True
# if Parameterized spatial resolution, set ParameterizeResolution to True
digiSIT.ParameterizeResolution = False
digiSIT.ParametersU = [2.29655e-03, 0.965899e-03, 0.584699e-03, 17.0856, 84.566, 12.4695e-03, -0.0643059, 0.168662, 1.87998e-03, 0.514452]
digiSIT.ParametersV = [1.44629e-02, 2.20108e-03, 1.03044e-02, 4.39195e+00, 3.29641e+00, 1.55167e+18, -5.41954e+01, 5.72986e+00, -6.80699e-03, 5.04095e-01]
#digiSIT.OutputLevel = DEBUG

digiSET = PlanarDigiAlg("SETDigi")
digiSET.IsStrip = False
digiSET.SimTrackHitCollection = "SETCollection"
digiSET.TrackerHitCollection = sethitname
digiSET.TrackerHitAssociationCollection = "SETTrackerHitAssociation"
digiSET.ResolutionU = [0.0072]
digiSET.ResolutionV = [0.086]
digiSET.UsePlanarTag = True
# if Parameterized spatial resolution, set ParameterizeResolution to True
digiSET.ParameterizeResolution = False
digiSET.ParametersU = [2.29655e-03, 0.965899e-03, 0.584699e-03, 17.0856, 84.566, 12.4695e-03, -0.0643059, 0.168662, 1.87998e-03, 0.514452]
digiSET.ParametersV = [1.44629e-02, 2.20108e-03, 1.03044e-02, 4.39195e+00, 3.29641e+00, 1.55167e+18, -5.41954e+01, 5.72986e+00, -6.80699e-03, 5.04095e-01]
#digiSET.OutputLevel = DEBUG

digiFTD = PlanarDigiAlg("FTDDigi")
digiFTD.IsStrip = False
digiFTD.SimTrackHitCollection = "FTDCollection"
digiFTD.TrackerHitCollection = ftdhitname
digiFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
digiFTD.ResolutionU = [0.003, 0.003, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072]
digiFTD.ResolutionV = [0.003, 0.003, 0.0072, 0.0072, 0.0072, 0.0072, 0.0072]
digiFTD.UsePlanarTag = True
# if Parameterized spatial resolution, set ParameterizeResolution to True
digiFTD.ParameterizeResolution = False
digiFTD.ParametersU = [2.29655e-03, 0.965899e-03, 0.584699e-03, 17.0856, 84.566, 12.4695e-03, -0.0643059, 0.168662, 1.87998e-03, 0.514452]
digiFTD.ParametersV = [1.44629e-02, 2.20108e-03, 1.03044e-02, 4.39195e+00, 3.29641e+00, 1.55167e+18, -5.41954e+01, 5.72986e+00, -6.80699e-03, 5.04095e-01]
#digiFTD.OutputLevel = DEBUG

from Configurables import DCHDigiAlg
digiDC = DCHDigiAlg("DCHDigi")
digiDC.DigiDCHitCollection = dchitname
# TODO: DCHDigiAlg need fixed, only WriteAna = True can run
digiDC.WriteAna = True
#digiDC.OutputLevel = DEBUG

# two strip tracker hits -> one space point
from Configurables import SpacePointBuilderAlg
spFTD = SpacePointBuilderAlg("FTDBuilder")
spFTD.TrackerHitCollection = ftdhitname
spFTD.TrackerHitAssociationCollection = "FTDTrackerHitAssociation"
spFTD.SpacePointCollection = ftdspname
spFTD.SpacePointAssociationCollection = "FTDSpacePointAssociation"
#spFTD.OutputLevel = DEBUG

# tracking
from Configurables import SiliconTrackingAlg
tracking = SiliconTrackingAlg("SiliconTracking")
tracking.HeaderCol = "EventHeader"
tracking.VTXHitCollection = vxdhitname
tracking.SITHitCollection = sithitname
tracking.FTDPixelHitCollection = ftdhitname
tracking.FTDSpacePointCollection = ftdspname
tracking.SITRawHitCollection = "NotNeedForPixelSIT"
tracking.FTDRawHitCollection = ftdhitname
tracking.UseSIT = True
tracking.SmoothOn = False
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
subset.RawTrackerHitCollections = [vxdhitname, sithitname, ftdhitname, ftdspname]
subset.TrackSubsetCollection = "SubsetTracks"
#subset.OutputLevel = DEBUG

#TODO: DC reconstruction, as preliminary, use Clupatra like as TPC
from Configurables import ClupatraAlg
clupatra = ClupatraAlg("Clupatra")
clupatra.TPCHitCollection = dchitname
clupatra.DistanceCut = 100.
clupatra.MaxDeltaChi2 = 100.
clupatra.Chi2Cut = 150.
#clupatra.OutputLevel = DEBUG

from Configurables import FullLDCTrackingAlg
full = FullLDCTrackingAlg("FullTracking")
full.VTXTrackerHits = vxdhitname
full.SITTrackerHits = sithitname
full.TPCTrackerHits = dchitname # add TPC or DC tracker hit here, if TPC or DC track is set by full.TPCTracks
full.SETTrackerHits = sethitname
full.FTDPixelTrackerHits = ftdhitname
full.FTDSpacePoints = ftdspname
full.SITRawHits     = "NotNeedForPixelSIT"
full.SETRawHits     = "NotNeedForPixelSET"
full.FTDRawHits     = ftdhitname
full.TPCTracks = "ClupatraTracks" # add standalone TPC or DC track here
full.SiTracks  = "SubsetTracks"
full.OutputTracks  = "MarlinTrkTracks"
full.SITHitToTrackDistance = 3.
full.SETHitToTrackDistance = 5.
#full.OutputLevel = DEBUG

#TODO: more reconstruction, PFA etc.

# output
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "CRD-o1-v01-SimRec00.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    #TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spFTD, digiDC, tracking, forward, subset, clupatra, full, out],
    TopAlg = [genalg, detsimalg, digiVXD, digiSIT, digiSET, digiFTD, spFTD, digiDC, tracking, forward, subset, full, out],
    EvtSel = 'NONE',
    EvtMax = 10,
    ExtSvc = [rndmengine, rndmgensvc, dsvc, evtseeder, geosvc, gearsvc, tracksystemsvc],
    HistogramPersistency = 'ROOT',
    OutputLevel = INFO
)
