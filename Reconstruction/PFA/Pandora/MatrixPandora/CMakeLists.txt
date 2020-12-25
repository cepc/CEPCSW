# Modules
gaudi_add_module(MatrixPandora
                 SOURCES src/PandoraMatrixAlg.cpp
                         src/MCParticleCreator.cpp
                         src/GeometryCreator.cpp
                         src/CaloHitCreator.cpp
                         src/TrackCreator.cpp
                         src/PfoCreator.cpp
                 LINK EventSeeder
                      GearSvc 
                      DetInterface
                      DataHelperLib
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

target_include_directories(MatrixPandora PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS MatrixPandora
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

