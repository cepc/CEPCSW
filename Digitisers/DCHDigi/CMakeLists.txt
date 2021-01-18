## Modules
gaudi_add_module(DCHDigi
                 SOURCES src/DCHDigiAlg.cpp
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
