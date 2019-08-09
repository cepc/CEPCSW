#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import DetSimSvc

detsimsvc = DetSimSvc("DetSimSvc")

from Configurables import DetSimAlg

detsimalg = DetSimAlg("DetSimAlg")

detsimalg.VisMacs = ["vis.mac"]

detsimalg.RunCmds = [
    "/tracking/verbose 1",
]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [detsimalg],
                EvtSel = 'NONE',
                EvtMax = 10,
)
