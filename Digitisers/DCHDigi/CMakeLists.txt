
find_package(k4FWCore)
find_package(CLHEP REQUIRED;CONFIG)
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED )

find_package(ROOT COMPONENTS MathCore Physics GenVector Geom REQUIRED)

set(DCHDigi_srcs
    src/DCHDigiAlg.cpp
)

## Modules
gaudi_add_module(DCHDigi
                 SOURCES ${DCHDigi_srcs}
                 LINK DetInterface
                      DetSegmentation
                      k4FWCore::k4FWCore 
                      Gaudi::GaudiKernel
                      Gaudi::GaudiAlgLib
                      ${CLHEP_LIBRARIES}
                      ${DD4hep_COMPONENT_LIBRARIES}
                      ${ROOT_LIBRARIES}
                      DetSegmentation
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
install(TARGETS DCHDigi
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
