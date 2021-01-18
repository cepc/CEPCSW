# Modules
gaudi_add_module(SimpleDigi 
                 SOURCES src/PlanarDigiAlg.cpp
                         src/TPCDigiAlg.cpp
                         src/voxel.cpp
                 LINK GearSvc
                      EventSeeder
                      TrackSystemSvcLib
                      DataHelperLib
                      k4FWCore::k4FWCore
                      Gaudi::GaudiKernel 
                      Gaudi::GaudiAlgLib
                      ${CLHEP_LIBRARIES}
                      ${GEAR_LIBRARIES}
                      ${GSL_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

install(TARGETS SimpleDigi
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
