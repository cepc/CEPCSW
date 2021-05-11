
# Modules
gaudi_add_module(GaudiPandora
                 SOURCES src/PandoraPFAlg.cpp
                         src/MCParticleCreator.cpp
                         src/GeometryCreator.cpp
                         src/CaloHitCreator.cpp
                         src/TrackCreator.cpp
                         src/PfoCreator.cpp
                         src/Utility.cpp
                 LINK EventSeeder
                      GearSvc
                      DataHelperLib
                      DetInterface
                      Gaudi::GaudiKernel
                      k4FWCore::k4FWCore
                      ${PandoraSDK_LIBRARIES}
                      ${LCContent_LIBRARIES}
                      ${CLHEP_LIBRARIES} 
                      ${ROOT_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      ${GEAR_LIBRARIES}
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

target_include_directories(GaudiPandora PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS GaudiPandora
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

