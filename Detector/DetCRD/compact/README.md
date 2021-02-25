# CRD detector models - Overview

The following CRD detector models are available in CEPCSW

| Model         |  Description                 | MainTracker |  Ecal   | Hcal | Status         |
| ------------- | -----------------------------|------------ |---------|------|----------------|
| CRD_o1_v01    | coil inside simulation model | DC          | crystal | RPC  | developing     |
| CRD_o1_v02    | pixel SIT                    | DC          | crystal | RPC  | developing     |
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
         - with Dirft Chamber + silicon layer between inner and outer chambers (DC + SIT34 + SET)
         - DC_outer_radius = 1716*mm
         - Detector/DetDriftChamber/src/driftchamber/DriftChamber.cpp
         - Detector/DetCEPCv4/src/tracker/SET_Simple_Planar_geo.cpp  
 - EndcapTracker
         - with silicon pestals (FTDPixel + FTDStrip)
         - Detector/DetCEPCv4/src/tracker/FTD_Simple_Staggered_geo.cpp
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
 - SIT:   strip -> pixel
 - SIT12->SIT1: r = 152.90  -> 140.00 mm (pixel SIT1)
 - SIT34->SIT2: between DCs -> 225.00 mm (pixel SIT2)
 - SET12: outside DC -> between DCs
 - SET34: outside DC (added)
 - compact files:
         - [./CRD_o1_v02/CRD_o1_v02.xml](./CRD_o1_v02/CRD_o1_v02.xml)
