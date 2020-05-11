
gaudi_subdir(DetSimCore v0r0)

gaudi_depends_on_subdirs(
    FWCore
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimCore_srcs
    src/DetSimAlg.cpp
    src/DetSimSvc.cpp
    src/DetectorConstruction.cpp
    src/PrimaryGeneratorAction.cpp
    src/G4PrimaryCnvTool.cpp
    src/ActionInitialization.cpp
    src/RunAction.cpp
    src/EventAction.cpp
    src/TrackingAction.cpp
    src/SteppingAction.cpp
)

message(" Geant4_LIBRARIES: ${Geant4_LIBRARIES}")

gaudi_add_module(DetSimCore ${DetSimCore_srcs}
    INCLUDE_DIRS
      # DetSimInterface
      # GaudiKernel
      # Geant4
      # plcio ${plcio_INCLUDE_DIRS} ${podio_INCLUDE_DIRS}
    LINK_LIBRARIES
      DetSimInterface
      GaudiKernel
      # Geant4
      # ${Geant4_LIBRARIES}
      DD4hep ${DD4hep_COMPONENT_LIBRARIES} 
      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
      # podio
)
