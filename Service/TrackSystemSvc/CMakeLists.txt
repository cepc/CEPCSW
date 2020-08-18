gaudi_subdir(TrackSystemSvc v0r0)

find_package(ROOT 6.14 REQUIRED COMPONENTS Matrix Physics)
find_package(CLHEP REQUIRED)
find_package(GEAR REQUIRED)
find_package(LCIO REQUIRED)
find_package(EDM4HEP REQUIRED)
#find_package(KalTest REQUIRED)
#find_package(KalDet REQUIRED)

gaudi_depends_on_subdirs(Service/GearSvc Detector/DetInterface Utilities/DataHelper Utilities/KalTest Utilities/KalDet)

set(TrackSystemSvc_srcs src/*.cpp)
set(TrackSystemSvcLib_srcs src/*.cc)

gaudi_install_headers(TrackSystemSvc)

#message( "${INCLUDE_DIRS}" )
#message( "${LINK_LIBRARIES}" )

gaudi_add_library(TrackSystemSvcLib ${TrackSystemSvcLib_srcs}
    PUBLIC_HEADERS TrackSystemSvc
    INCLUDE_DIRS GaudiKernel ROOT CLHEP gear ${LCIO_INCLUDE_DIRS} ${EDM4HEP_INCLUDE_DIRS}
    LINK_LIBRARIES DataHelperLib KalTestLib KalDetLib GaudiKernel ROOT CLHEP $ENV{GEAR}/lib/libgear.so $ENV{GEAR}/lib/libgearsurf.so ${LCIO_LIBRARIES}
     -Wl,--no-as-needed
     EDM4HEP::edm4hep EDM4HEP::edm4hepDict
     -Wl,--as-needed
)

gaudi_add_module(TrackSystemSvc ${TrackSystemSvc_srcs}
				INCLUDE_DIRS GaudiKernel gear 
				LINK_LIBRARIES TrackSystemSvcLib GaudiKernel $ENV{GEAR}/lib/libgear.so
)
