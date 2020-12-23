# gaudi_subdir(TrackSystemSvc v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(ROOT 6.14 REQUIRED COMPONENTS Matrix Physics)
find_package(GEAR REQUIRED)
find_package(LCIO REQUIRED)
find_package(EDM4HEP REQUIRED)
#find_package(KalTest REQUIRED)
#find_package(KalDet REQUIRED)
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)

# gaudi_depends_on_subdirs(Service/GearSvc Detector/DetInterface Utilities/DataHelper Utilities/KalTest Utilities/KalDet)

set(TrackSystemSvc_srcs src/*.cpp)
set(TrackSystemSvcLib_srcs src/*.cc)

# gaudi_install_headers(TrackSystemSvc)

#message( "${INCLUDE_DIRS}" )
#message( "${LINK_LIBRARIES}" )

gaudi_add_library(TrackSystemSvcLib
                  SOURCES ${TrackSystemSvcLib_srcs}
                  LINK DataHelperLib 
                       KalTestLib 
                       KalDetLib 
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
