
# Modules
gaudi_add_module(TrackInspect
                 SOURCES src/TrackInspectAlg.cpp
                 LINK DataHelperLib 
                      Gaudi::GaudiKernel 
                      EDM4HEP::edm4hep 
                      k4FWCore::k4FWCore
                      ${ROOT_LIBRARIES}
                      ${CLHEP_LIBRARIES}
)

install(TARGETS TrackInspect
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
