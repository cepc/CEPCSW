#!/usr/bin/env python

from Gaudi.Configuration import *

from Configurables import HelloAlg

helloalg = HelloAlg("helloAlg")
helloalg.MyInt = 42

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [helloalg],
                EvtSel = 'NONE',
                EvtMax = 10,
)
