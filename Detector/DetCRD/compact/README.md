# CRD detector models - Overview

The following CRD detector models are available in CEPCSW

| Model         |  Description                 | MainTracker |  Ecal   | Hcal | Status         |
| ------------- | -----------------------------|------------ |---------|------|----------------|
| CRD_i1_v01    | coil inside simulation model | DC          | crystal | RPC  | developing     |
| ------------- | -----------------------------|-------------|---------|------|----------------|

## Details

### CRD_i1_v01
 - coil inside CRD model
 - BeamPipe
         - with center pipe + crotch link to doubly-pipe
 - Vertex
         - with silicon ladders (VXD + SIT)
 - MainTracker
         - with Dirft Chamber + silicon layer between inner and outer chambers
         - DC_outer_radius = 1800*mm
 - EndcapTracker
         - with silicon pestals (FTDPixel + FTDStrip)
 - Ecal
         - with crystal 
 - Hcal
         - with scintillator **and** RPC readout
         - creates two sets of hit collections
 - compact files:
         - [./CRD_i1_v01/CRD_i1_v01.xml](./CRD_i1_v01/CRD_i1_v01.xml)


