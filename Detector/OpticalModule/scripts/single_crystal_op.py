#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import K4DataSvc
dsvc = K4DataSvc("EventDataSvc")

from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
# rndmengine = HepRndm__Engine_CLHEP__RanluxEngine_() # The default engine in Gaudi
rndmengine = HepRndm__Engine_CLHEP__HepJamesRandom_() # The default engine in Geant4
rndmengine.SetSingleton = True
rndmengine.Seeds = [42]


geometry_option = "SingleCrystalOptical.xml"

if not os.getenv("OPTICALMODULEROOT"):
    print("Can't find the geometry. Please setup envvar OPTICALMODULEROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("OPTICALMODULEROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeoSvc
geosvc = GeoSvc("GeoSvc")
geosvc.compact = geometry_path

from Configurables import GenAlgo
from Configurables import GtGunTool
from Configurables import StdHepRdr
from Configurables import SLCIORdr
from Configurables import HepMCRdr
from Configurables import GenPrinter

gun = GtGunTool("GtGunTool")
#gun.Particles = ["opticalphoton"]
gun.Particles = ["mu-"]
gun.EnergyMins = [10]#2.6e-9] # GeV
gun.EnergyMaxs = [10]#2.6e-9]
gun.ThetaMins = [90] # rad; 45deg
gun.ThetaMaxs = [90] # rad; 45deg
gun.PhiMins = [0] # rad; 0deg
gun.PhiMaxs = [0] # rad; 360deg

genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtGunTool"]

##############################################################################
# Detector Simulation
##############################################################################
from Configurables import DetSimSvc
detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg
detsimalg = DetSimAlg("DetSimAlg")
# detsimalg.VisMacs = ["vis.mac"]
detsimalg.RunCmds = [
    "/physics_lists/factory/addOptical",
    "/process/optical/verbose 0"
]
detsimalg.PhysicsList = "FTFP_BERT"
detsimalg.AnaElems = ["Edm4hepWriterAnaElemTool"]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import PodioOutput
out = PodioOutput("Write")
out.filename = "test.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [genalg, detsimalg, out],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [rndmengine, dsvc, geosvc],
    HistogramPersistency = "ROOT",
    OutputLevel=INFO
)
