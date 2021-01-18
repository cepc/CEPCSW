# gaudi_subdir(TrackSystemSvc v0r0)

gaudi_add_library(TrackSystemSvcLib
                  SOURCES src/HelixFit.cc
                          src/HelixTrack.cc
                          src/IMarlinTrack.cc
                          src/IMarlinTrkSystem.cc
                          src/LCIOTrackPropagators.cc
                          src/MarlinKalTest.cc
                          src/MarlinKalTestTrack.cc
                          src/MarlinTrkUtils.cc

                  LINK DataHelperLib 
                       KalTestLib 
                       KalDetLib 
                       GearSvc
                       Gaudi::GaudiKernel 
                       ${ROOT_LIBRARIES}
                       ${CLHEP_LIBRARIES} 
                       ${GEAR_LIBRARIES}
                       ${LCIO_LIBRARIES}
                       EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

gaudi_add_module(TrackSystemSvc 
                 SOURCES src/TrackSystemSvc.cpp
		 LINK TrackSystemSvcLib 
                      Gaudi::GaudiKernel
                      ${GEAR_LIBRARIES}
)

install(TARGETS TrackSystemSvcLib TrackSystemSvc
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

