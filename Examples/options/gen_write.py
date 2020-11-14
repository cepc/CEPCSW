#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc")

from Configurables import GenAlgo

read = GenAlgo("read")
#######################################
#support format: stdhep, slcio, hepmc #
#######################################
#read.Input = "/junofs/users/wxfang/CEPC/CEPCOFF/TestExample/stdhep/nnh_e2e2.e0.p0.00001.stdhep"
#read.FileFormat = "stdhep"
read.Input = "/junofs/users/wxfang/CEPC/whizard_apply/ee/ee.slcio"
read.FileFormat = "slcio"
read.PrintEvent = True # true for printing mc info
read.WriteFile = True  # true for writting info to root

from Configurables import PodioOutput
out = PodioOutput("out")
out.filename = "test.root" #name of output root file
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [read, out],
                EvtSel = 'NONE',
                EvtMax = 11,
                ExtSvc=[dsvc],
                OutputLevel=DEBUG
)
