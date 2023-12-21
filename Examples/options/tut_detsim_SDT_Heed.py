#!/usr/bin/env python

import os
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
#geometry_option = "CepC_v4.xml"
det_root = "DETDRIFTCHAMBERROOT"
#det_root = "DETCEPCV4ROOT"#"DETDRIFTCHAMBERROOT"
if not os.getenv(det_root):
    print("Can't find the geometry. Please setup envvar %s."%det_root )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv(det_root), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
print('geometry_path=',geometry_path)
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


gun.EnergyMins = [10] # GeV
gun.EnergyMaxs = [10] # GeV

gun.ThetaMins = [80] # rad; 45deg
gun.ThetaMaxs = [90.] # rad; 45deg

gun.PhiMins = [0] # rad; 0deg
gun.PhiMaxs = [360.] # rad; 360deg

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

from Configurables import DummyFastSimG4Tool
dummy_fastsim_tool = DummyFastSimG4Tool("DummyFastSimG4Tool")

detsimalg.FastSimG4Tools = [
#    "DummyFastSimG4Tool"
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

#dedxoption = "DummyDedxSimTool"
#dedxoption = "BetheBlochEquationDedxSimTool"
dedxoption = "TrackHeedSimTool"

driftchamber_sensdettool.DedxSimTool = dedxoption

from Configurables import DummyDedxSimTool
from Configurables import BetheBlochEquationDedxSimTool
from Configurables import TrackHeedSimTool

if dedxoption == "DummyDedxSimTool":
    dedx_simtool = DummyDedxSimTool("DummyDedxSimTool")
elif dedxoption == "BetheBlochEquationDedxSimTool":
    dedx_simtool = BetheBlochEquationDedxSimTool("BetheBlochEquationDedxSimTool")
    dedx_simtool.material_Z = 2
    dedx_simtool.material_A = 4
    dedx_simtool.scale = 10
    dedx_simtool.resolution = 0.0001
elif dedxoption == "TrackHeedSimTool":
    dedx_simtool = TrackHeedSimTool("TrackHeedSimTool")
    dedx_simtool.only_primary = False#True
    dedx_simtool.use_max_step = True#True
    dedx_simtool.max_step = 1#mm
    #dedx_simtool.he   = 50
    #dedx_simtool.isob = 50
    #dedx_simtool.gas_file ="/junofs/users/wxfang/MyGit/tmp/check_G4FastSim_20210121/CEPCSW/Digitisers/DigiGarfield/He_50_isobutane_50.gas" 
    dedx_simtool.he   = 90
    dedx_simtool.isob = 10
    #dedx_simtool.gas_file ="/junofs/users/wxfang/MyGit/tmp/check_G4FastSim_20210121/CEPCSW/Digitisers/DigiGarfield/he_90_isobutane_10.gas" 
    #dedx_simtool.IonMobility_file ="/junofs/users/wxfang/MyGit/tmp/check_G4FastSim_20210121/CEPCSW/Digitisers/DigiGarfield/IonMobility_He+_He.txt" 
    dedx_simtool.gas_file         ="he_90_isobutane_10.gas"
    dedx_simtool.IonMobility_file ="IonMobility_He+_He.txt"
    dedx_simtool.save_mc = True
    dedx_simtool.debug = False
    dedx_simtool.sim_pulse = True
    #dedx_simtool.model='/junofs/users/wxfang/MyGit/tmp/fork_cepcsw_20220418/CEPCSW/Digitisers/SimCurrentONNX/src/model_90He10C4H10_18mm.onnx'
    dedx_simtool.model='model_90He10C4H10_18mm.onnx'
    dedx_simtool.batchsize = 100

##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "detsim_heed.root"
out.outputCommands = ["keep *"]

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg, out],
                EvtSel = 'NONE',
                EvtMax = 20,
                ExtSvc = [rndmengine, rndmgensvc, dsvc, geosvc],
                OutputLevel=INFO
)
