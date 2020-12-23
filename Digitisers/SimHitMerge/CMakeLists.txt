
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore REQUIRED)

# Modules
gaudi_add_module(SimHitMerge
                 SOURCES src/SimHitMergeAlg.cpp
                 LINK DetInterface
                      k4FWCore::k4FWCore
                      Gaudi::GaudiKernel
                      Gaudi::GaudiAlgLib
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

install(TARGETS SimHitMerge
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

