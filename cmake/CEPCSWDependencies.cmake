#[[

Find all the dependencies here, so in each package user don't need to find the packages again.

- CLHEP
- DD4hep
- EDM4hep
- Gaudi
- Geant4
- GEAR
- GSL
- HepMC
- k4FWCore
- LCContent
- LCIO
- PandoraSDK
- podio
- ROOT
#]]

find_package(CLHEP REQUIRED;CONFIG)
find_package(DD4hep COMPONENTS DDCore DDG4 DDParsers DDRec REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED)
find_package(HepMC)
find_package(k4FWCore REQUIRED)
find_package(LCContent REQUIRED)
find_package(LCIO REQUIRED)
find_package(PandoraSDK REQUIRED)
find_package(podio REQUIRED)
find_package(ROOT COMPONENTS EG Graf Graf3d Gpad MathCore Net RIO Tree TreePlayer REQUIRED)
find_package(GenFit)
