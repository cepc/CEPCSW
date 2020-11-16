#!/usr/bin/env python

from Gaudi.Configuration import *
from Configurables import K4DataSvc
#dsvc = K4DataSvc("EventDataSvc", input="detsim_ECAL_gamma_10000evt.root")
dsvc = K4DataSvc("EventDataSvc", input="ECALonly_gamma_30degree_10000evt.root")

#########################################################################
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
########################################################################

from Configurables import PodioInput
podioinput = PodioInput("PodioReader", collections=[
    "MCParticle",
    "EcalBarrelCollection",
    "EcalEndcapsCollection"
])


############################################################
from Configurables import GearSvc
gearSvc  = GearSvc("GearSvc")
gearSvc.GearXMLFile = "Detector/DetCEPCv4/compact/FullDetGear.xml"

############################################################
from Configurables import SimHitMergeAlg
simHitMerge = SimHitMergeAlg("SimHitMergeAlg")
simHitMerge.InputCollections=["EcalBarrelCollection", "EcalEndcapsCollection"]
simHitMerge.OutputCollections=["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged"]
############################################################

from Configurables import G2CDArborAlg
caloDigi = G2CDArborAlg("G2CDArborAlg")
caloDigi.ReadLCIO = False 
#caloDigi.CalibrECAL = [48.16, 96.32]
caloDigi.CalibrECAL = [46.538, 93.0769]
caloDigi.ECALCollections = ["EcalBarrelCollectionMerged", "EcalEndcapsCollectionMerged"]
caloDigi.ECALReadOutNames= ["EcalBarrelCollection", "EcalEndcapsCollection"]
caloDigi.DigiECALCollection = ["EcalBarrel", "EcalEndcaps"]
caloDigi.HCALCollections = []
caloDigi.HCALReadOutNames= []
caloDigi.DigiHCALCollection = []
caloDigi.EventReportEvery = 100
##############################################################################


#from Configurables import DumpDigiIDAlgv1
#dumpIDAlg = DumpDigiIDAlgv1("DumpDigiIDAlgv1")
#dumpIDAlg.ColName = "EcalEndcaps"
#dumpIDAlg.Readout = "EcalEndcapsCollection"
##dumpIDAlg.Output = "id_ECAL_gamma_digi_10000evt_newCalibv2_fromMergedSimHit.root"
#dumpIDAlg.Output = "id_ECALonly_30degree.root"
##############################################################################

from Configurables import PodioOutput
out = PodioOutput("outputalg")
out.filename = "digi_fix5GeV_EcalOnlyBF.root"
out.outputCommands = ["keep *"]

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(
TopAlg = [podioinput, simHitMerge, caloDigi, out],
    EvtSel = 'NONE',
    EvtMax = -1,
    ExtSvc = [dsvc, gearSvc],
    OutputLevel=INFO
)
