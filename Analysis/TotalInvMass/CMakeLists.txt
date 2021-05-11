
gaudi_add_module(TotalInvMass
                 SOURCES src/TotalInvMass.cc
                 LINK k4FWCore::k4FWCore 
                      Gaudi::GaudiKernel 
                      Gaudi::GaudiAlgLib 
                      ${CLHEP_LIBRARIES} 
                      ${GSL_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ${ROOT_LIBRARIES}
)

install(TARGETS TotalInvMass
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
