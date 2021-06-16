# CRD detector models - Overview

The following CRD detector models are available in CEPCSW

| Model         |  Description                 | MainTracker |  Ecal   | Hcal | Status         |
| ------------- | -----------------------------|------------ |---------|------|----------------|
| CRD_o1_v01    | coil inside simulation model | SIT+DC+SET  | crystal | RPC  | developing     |
| CRD_o1_v02    | pixel SET                    | SIT+DC+SET  | crystal | RPC  | developing     |
| ------------- | -----------------------------|-------------|---------|------|----------------|
 
## Details

### CRD_o1_v01 (to update)
 - coil inside CRD model
 - BeamPipe
         - with center pipe + crotch link to doubly-pipe
         - Detector/DetCRD/src/Other/CRDBeamPipe_v01_geo.cpp
 - Vertex
         - with silicon ladders (VXD + SIT12)
         - Detector/DetCEPCv4/src/tracker/VXD04_geo.cpp
         - Detector/DetCEPCv4/src/tracker/SIT_Simple_Planar_geo.cpp
 - MainTracker
         - with Dirft Chamber + pixel silicon layer between inner and outer chambers (SIT1234 + DC + SET)
         - DC_outer_radius = 1800*mm for maximum sensitive gas 
         - Detector/DetDriftChamber/src/driftchamber/DriftChamber.cpp
         - Detector/DetCEPCv4/src/tracker/SIT_Simple_Pixel_geo.cpp  
 - EndcapTracker
         - with silicon pestals (FTDPixel)
         - Detector/DetCRD/src/Tracker/SiTrackerSkewRing_v01_geo.cpp
 - Ecal
         - with crystal 
         - Detector/DetCRD/src/Calorimeter/CRDEcal.cpp
         - Endcap (TODO)
 - Hcal
         - with RPC readout
         - creates two sets of hit collections
 - Coil
         - CEPC_v4 like
 - Yoke
         - CEPC_v4 like
 - compact files:
         - [./CRD_o1_v01/CRD_o1_v01.xml](./CRD_o1_v01/CRD_o1_v01.xml)

### CRD_o1_v02 (to update)
 - based on CRD_o1_v01
 - strip SET: double layers
 - compact files:
         - [./CRD_o1_v02/CRD_o1_v02.xml](./CRD_o1_v02/CRD_o1_v02.xml)
