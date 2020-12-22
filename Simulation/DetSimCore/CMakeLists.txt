
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

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

gaudi_add_module(DetSimCore 
                 SOURCES ${DetSimCore_srcs}
                 LINK DetSimInterface
                      Gaudi::GaudiKernel
                      ${Geant4_LIBRARIES}
                      ${DD4hep_COMPONENT_LIBRARIES} 
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

install(TARGETS DetSimCore
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


