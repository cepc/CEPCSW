#!/usr/bin/env python

import os
print(os.environ["DD4HEP_LIBRARY_PATH"])
import sys
# sys.exit(0)

from Gaudi.Configuration import *

##############################################################################
# Event Data Svc
##############################################################################
from Configurables import CEPCDataSvc
dsvc = CEPCDataSvc("EventDataSvc")

##############################################################################
# Physics Generator
##############################################################################
from Configurables import GenAlgo

genalg = GenAlgo("read")
genalg.Input = "/cefs/data/stdhep/CEPC250/2fermions/E250.Pbhabha.e0.p0.whizard195/bhabha.e0.p0.00001.stdhep"
genalg.FileFormat = "stdhep"
genalg.PrintEvent = False # true for printing mc info

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

# geometry_option = "CepC_v4-onlyTracker.xml"
geometry_option = "CepC_v4-onlyVXD.xml"

if os.getenv("DETCEPCV4ROOT"):
    example_dettool.detxml = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
else:
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

##############################################################################
# POD I/O
##############################################################################
from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "test-detsim10.root"
out.outputCommands = ["keep *"]

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg, out],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
)
