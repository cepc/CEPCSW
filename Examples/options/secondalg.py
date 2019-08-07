#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import FirstSvc

firstsvc = FirstSvc("FirstSvc")

from Configurables import SecondAlg

secondalg = SecondAlg("secondAlg")

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [secondalg],
                EvtSel = 'NONE',
                EvtMax = 10,
)
