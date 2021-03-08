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

# geometry_option = "CepC_v4-onlyTracker.xml"
geometry_option = "det.xml"

if not os.getenv("DETDRIFTCHAMBERROOT"):
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETDRIFTCHAMBERROOT"), "compact", geometry_option)
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
# gun.Particles = ["pi+"]
# gun.EnergyMins = [100.] # GeV
# gun.EnergyMaxs = [100.] # GeV

gun.Particles = ["e-"]

# gun.PositionXs = [100.] # mm
# gun.PositionYs = [100.] # mm
# gun.PositionZs = [0.] # mm


gun.EnergyMins = [10.] # GeV
gun.EnergyMaxs = [10.] # GeV

gun.ThetaMins = [90] # rad; 45deg
gun.ThetaMaxs = [90] # rad; 45deg

gun.PhiMins = [90] # rad; 0deg
gun.PhiMaxs = [90] # rad; 360deg
#gun.PhiMins = [0] # rad; 0deg
#gun.PhiMaxs = [360] # rad; 360deg

# stdheprdr = StdHepRdr("StdHepRdr")
# stdheprdr.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"

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

# from Configurables import ExampleAnaElemTool
# example_anatool = ExampleAnaElemTool("ExampleAnaElemTool")

from Configurables import DetSimAlg

detsimalg = DetSimAlg("DetSimAlg")
detsimalg.RandomSeeds = seed

if int(os.environ.get("VIS", 0)):
    detsimalg.VisMacs = ["vis.mac"]

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
from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
ntsvc.Output = ["MyTuples DATAFILE='DCH_digi_ana.root' OPT='NEW' TYP='ROOT'"]

##############################################################################
# DCHDigiAlg
##############################################################################
from Configurables import DCHDigiAlg
dCHDigiAlg = DCHDigiAlg("DCHDigiAlg")
dCHDigiAlg.readout = "DriftChamberHitsCollection"
dCHDigiAlg.drift_velocity = 40#um/ns
dCHDigiAlg.mom_threshold = 0 #GeV
dCHDigiAlg.SimDCHitCollection = "DriftChamberHitsCollection"
dCHDigiAlg.DigiDCHitCollection = "DigiDCHitsCollection"
dCHDigiAlg.AssociationCollection = "DCHAssociationCollectio"
dCHDigiAlg.WriteAna  = True

##############################################################################
# TruthTrackerAlg
##############################################################################
from Configurables import TruthTrackerAlg
truthTrackerAlg = TruthTrackerAlg("TruthTrackerAlg")
truthTrackerAlg.DCHitAssociationCollection="DCHAssociationCollectio"
# truthTrackerAlg.debug = 1

##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "truthRec_DCH.root"
out.outputCommands = ["keep *"]

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg, dCHDigiAlg, truthTrackerAlg, out],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [rndmengine, rndmgensvc, dsvc, geosvc],
                HistogramPersistency = "ROOT",
                OutputLevel=INFO
)
