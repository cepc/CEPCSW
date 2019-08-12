#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import DetSimSvc

detsimsvc = DetSimSvc("DetSimSvc")

# from Configurables import ExampleAnaElemTool
# example_anatool = ExampleAnaElemTool("ExampleAnaElemTool")

from Configurables import DetSimAlg

detsimalg = DetSimAlg("DetSimAlg")

detsimalg.VisMacs = ["vis.mac"]

detsimalg.RunCmds = [
    "/tracking/verbose 1",
]
detsimalg.AnaElems = [
    # example_anatool.name()
    "ExampleAnaElemTool"
]
detsimalg.RootDetElem = "WorldDetElemTool"

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [detsimalg],
                EvtSel = 'NONE',
                EvtMax = 10,
)
