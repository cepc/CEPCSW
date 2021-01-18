
# Modules
gaudi_add_module(DCHDedx
                 SOURCES src/RecDCHDedxAlg.cpp
                 LINK DetSimInterface
                      DetInterface
                      DetSegmentation
                      k4FWCore::k4FWCore
                      Gaudi::GaudiAlgLib 
                      Gaudi::GaudiKernel 
                      DetSegmentation
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

install(TARGETS DCHDedx
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
