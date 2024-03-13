#!/usr/bin/env python
#Author: Zhan Li <lizhan@ihep.ac.cn>
#Created [2024-03-07 Thu 14:53]

import os
import sys

import Gaudi.Configuration
from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_, k4DataSvc, GeomSvc
from Configurables import TimeProjectionChamberSensDetTool
from Configurables import GenAlgo
from Configurables import GtGunTool
from Configurables import StdHepRdr
from Configurables import SLCIORdr
from Configurables import HepMCRdr
from Configurables import GenPrinter
from Configurables import GtBeamBackgroundTool
from Configurables import DetSimSvc
from Configurables import DetSimAlg
from Configurables import AnExampleDetElemTool
from Configurables import PodioOutput
from Configurables import ApplicationMgr

seed = [42]

rndmengine = Gaudi.Configuration.HepRndm__Engine_CLHEP__HepJamesRandom_("RndmGenSvc.Engine")
rndmengine.SetSingleton = True
rndmengine.Seeds = seed

rndmgensvc = RndmGenSvc("RndmGenSvc")
rndmgensvc.Engine = rndmengine.name()


dsvc = k4DataSvc("EventDataSvc")
#geometry_option = "CepC_v4-onlyVXD.xml"
geometry_option = "CepC_v4_onlyTracker.xml"
#geometry_option = "CepC_v4.xml"

geometry_path = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path

#Previously I do not have these 2 lines
tpc_sensdettool = TimeProjectionChamberSensDetTool("TimeProjectionChamberSensDetTool")
tpc_sensdettool.TypeOption = 1


# Physics Generator
bg = GtBeamBackgroundTool("GtBeamBackgroundTool")
bg.InputFileMap = {"default":"/scratchfs/atlas/lizhan/cepc/CEPCSW/ToCEPCSWsingle.out"}
#bg.InputFileMap = {"default":"/cefs/higgs/shihy/tools/CEPCSW/CEPCSW/Test/0/ToCEPCSW-1.out"}
bg.InputBeamEnergyMap = {"default":120}
bg.RotationAlongYMap = {"default":16.5e-3}


genprinter = GenPrinter("GenPrinter")

genalg = GenAlgo("GenAlgo")
genalg.GenTools = ["GtBeamBackgroundTool"]

detsimsvc = DetSimSvc("DetSimSvc")

detsimalg = DetSimAlg("DetSimAlg")
detsimalg.RandomSeeds = seed


detsimalg.RunCmds = []
detsimalg.AnaElems = [
    "Edm4hepWriterAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"

example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")


# POD I/O
out = PodioOutput("outputalg")
out.filename = "test-SIT.root"
out.outputCommands = ["keep *"]

ApplicationMgr( TopAlg = [genalg, detsimalg, out],
                EvtSel = 'NONE',
                EvtMax = 1,
                ExtSvc = [rndmengine, rndmgensvc, dsvc, geosvc],
)