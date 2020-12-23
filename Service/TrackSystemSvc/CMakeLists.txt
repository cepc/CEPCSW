# gaudi_subdir(TrackSystemSvc v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(ROOT 6.14 REQUIRED COMPONENTS Matrix Physics)
find_package(GEAR REQUIRED)
find_package(LCIO REQUIRED)
find_package(EDM4HEP REQUIRED)
#find_package(KalTest REQUIRED)
#find_package(KalDet REQUIRED)
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)

set(TrackSystemSvc_srcs src/TrackSystemSvc.cpp)
set(TrackSystemSvcLib_srcs 
      src/HelixFit.cc
      src/HelixTrack.cc
      src/IMarlinTrack.cc
      src/IMarlinTrkSystem.cc
      src/LCIOTrackPropagators.cc
      src/MarlinKalTest.cc
      src/MarlinKalTestTrack.cc
      src/MarlinTrkUtils.cc
)

# gaudi_install_headers(TrackSystemSvc)

#message( "${INCLUDE_DIRS}" )
#message( "${LINK_LIBRARIES}" )

gaudi_add_library(TrackSystemSvcLib
                  SOURCES ${TrackSystemSvcLib_srcs}
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
                 SOURCES ${TrackSystemSvc_srcs}
		 LINK TrackSystemSvcLib 
                      Gaudi::GaudiKernel
                      ${GEAR_LIBRARIES}
)

install(TARGETS TrackSystemSvcLib TrackSystemSvc
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

