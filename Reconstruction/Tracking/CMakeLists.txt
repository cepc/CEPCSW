# Modules
gaudi_add_module(Tracking
                 SOURCES src/Clupatra/ClupatraAlg.cpp
                         src/Clupatra/clupatra_new.cpp
                         src/FullLDCTracking/FullLDCTrackingAlg.cpp
                         src/TruthTracker/TruthTrackerAlg.cpp
                 LINK GearSvc
                      EventSeeder
                      TrackSystemSvcLib
                      DetSegmentation
                      Gaudi::GaudiAlgLib
                      Gaudi::GaudiKernel
                      ${GEAR_LIBRARIES}
                      ${GSL_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      DetSegmentation
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

target_include_directories(Tracking PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS Tracking
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
