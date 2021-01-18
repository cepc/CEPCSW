# Modules
gaudi_add_module(G2CDArbor
                 SOURCES src/G2CDArborAlg.cpp
                 LINK k4FWCore::k4FWCore
                      GearSvc
                      DetInterface
                      Gaudi::GaudiKernel
                      Gaudi::GaudiAlgLib 
                      ${CLHEP_LIBRARIES}
                      ${GEAR_LIBRARIES} 
                      ${GSL_LIBRARIES} 
                      ${LCIO_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
install(TARGETS G2CDArbor
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

