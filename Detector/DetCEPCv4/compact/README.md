# CEPC_v4 detector models - Overview
# inplement according database models03 (MokkaC), by using lcgeo detector construction 

The following CEPC_v4 detector models are available in CEPCSW

| Model         |  Description                 | MainTracker |  Ecal   | Hcal | Status         |
| ------------- | -----------------------------|------------ |---------|------|----------------|
| CEPC_v4       | following MokkaC's version   | TPC         | Si-W    | RPC  | implementing   |
| ------------- | -----------------------------|-------------|---------|------|----------------|

## Details

### CEPC_v4
 - BeamPipe
         - Z from 0 to 700mm same as MokkaC's CEPC_v4
 - Vertex
         - with silicon ladders (VXD + SIT)
 - MainTracker
         - TPC
         - TPC_outer_radius = 1808*mm
 - EndcapTracker
         - with silicon pestals (FTDPixel + FTDStrip)
 - Ecal (un-implemented)
         - with si-W calorimeter
 - Hcal (un-implemented)
         - with scintillator **and** RPC readout
         - creates two sets of hit collections
 - compact files:
         - only Tracker [./CEPC_v4.xml](./CEPC_v4.xml)
