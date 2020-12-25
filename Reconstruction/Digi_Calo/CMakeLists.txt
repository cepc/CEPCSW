
## Modules
gaudi_add_module(Digi_Calo
                 SOURCES src/CaloDigiAlg.cpp
                 LINK DetInterface
                      k4FWCore::k4FWCore
                      Gaudi::GaudiKernel
                      Gaudi::GaudiAlgLib
                      ${CLHEP_LIBRARIES} 
                      ${DD4hep_COMPONENT_LIBRARIES} 
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
install(TARGETS Digi_Calo
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

