#!/usr/bin/env python

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
genalg.Input = "/junofs/users/wxfang/CEPC/whizard_apply/ee/ee.slcio"
genalg.FileFormat = "slcio"
genalg.PrintEvent = True # true for printing mc info
genalg.WriteFile = True  # true for writting info to root

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
    "/tracking/verbose 1",
]
detsimalg.AnaElems = [
    # example_anatool.name()
    "ExampleAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"

from Configurables import AnExampleDetElemTool
example_dettool = AnExampleDetElemTool("AnExampleDetElemTool")
example_dettool.detxml = "/cvmfs/sft.cern.ch/lcg/releases/DD4hep/01-08-c926f/x86_64-slc6-gcc62-opt/DDDetectors/compact/SiD.xml"

##############################################################################
# ApplicationMgr
##############################################################################

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [genalg, detsimalg],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc],
)
