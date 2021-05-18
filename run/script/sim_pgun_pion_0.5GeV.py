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
from Configurables import CEPCDataSvc
dsvc = CEPCDataSvc("EventDataSvc")


##############################################################################
# Geometry Svc
##############################################################################

geometry_option = "CepC_v4-onlyTracker.xml"
# geometry_option = "CepC_v4-onlyVXD.xml"
#geometry_option = "CepC_v4-onlyVXD_DCH.xml"

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
gun.Particles = ["pi+"]
# gun.Particles = ["mu+"]
#gun.Particles = ["proton"]
# gun.Energies = [100.] # GeV
# gun.Energies = [10.0] # GeV
gun.Energies = [0.5] # GeV

# gun.ThetaMins = [] # rad; 45deg
# gun.ThetaMaxs = [] # rad; 45deg

gun.ThetaMins = [ 85.0 ] # rad; 45deg
gun.ThetaMaxs = [ 85.0 ] # rad; 45deg


gun.PhiMins = [] # rad; 0deg
gun.PhiMaxs = [] # rad; 360deg

stdheprdr = StdHepRdr("StdHepRdr")
#stdheprdr.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"
stdheprdr.Input = "/cefs/data/stdhep/CEPC240/2fermions/E240.Pe2e2.e0.p0.whizard195/e2e2.e0.p0.00001.stdhep"

# lciordr = SLCIORdr("SLCIORdr")
# lciordr.Input = "/cefs/data/stdhep/lcio250/signal/Higgs/E250.Pbbh.whizard195/E250.Pbbh_X.e0.p0.whizard195/Pbbh_X.e0.p0.00001.slcio"

# hepmcrdr = HepMCRdr("HepMCRdr")
# hepmcrdr.Input = "example_UsingIterators.txt"

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtGunTool"]
##genalg.GenTools = ["StdHepRdr"]
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
    "ExampleAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")


##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
#out.filename = "test-detsim10.root"
out.filename = "../run/output/pion_0.5GeV_85deg_2000event.root"
out.outputCommands = ["keep *"]

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg, out],
                EvtSel = 'NONE',
                EvtMax = 2000,
                ExtSvc = [rndmengine, dsvc, geosvc],
)
