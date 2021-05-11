#!/usr/bin/env python

from Gaudi.Configuration import *

##############################################################################
# Geometry Svc
##############################################################################

# geometry_option = "CepC_v4-onlyTracker.xml"
# geometry_option = "CepC_v4-onlyVXD.xml"
geometry_option = "CepC_v4-onlyECAL.xml"

if not os.getenv("DETCEPCV4ROOT"):
    print("Can't find the geometry. Please setup envvar DETCEPCV4ROOT." )
    sys.exit(-1)

geometry_path = os.path.join(os.getenv("DETCEPCV4ROOT"), "compact", geometry_option)
if not os.path.exists(geometry_path):
    print("Can't find the compact geometry file: %s"%geometry_path)
    sys.exit(-1)

from Configurables import GeomSvc
geosvc = GeomSvc("GeomSvc")
geosvc.compact = geometry_path

##############################################################################
# Event Data Svc
##############################################################################

from Configurables import k4DataSvc
dsvc = k4DataSvc("EventDataSvc", input="test-detsim10.root")

##############################################################################
# NTuple Svc
##############################################################################

from Configurables import NTupleSvc
ntsvc = NTupleSvc("NTupleSvc")
ntsvc.Output = ["MyTuples DATAFILE='result.root' OPT='NEW' TYP='ROOT'"]

##############################################################################
# DumpAlg
##############################################################################

from Configurables import DumpIDAlg
alg = DumpIDAlg("DumpAlg")

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
        "EcalBarrelCollection"
    ])

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [podioinput, alg],
                EvtSel = 'NONE',
                EvtMax = 10,
                ExtSvc = [dsvc, ntsvc],
                HistogramPersistency = "ROOT",
                OutputLevel=DEBUG
)
