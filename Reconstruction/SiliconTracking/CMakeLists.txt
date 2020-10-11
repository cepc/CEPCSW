gaudi_subdir(SiliconTracking v0r0)

find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED )
find_package(EDM4HEP REQUIRED ) 

gaudi_depends_on_subdirs(
    Service/GearSvc
    Service/EventSeeder
    Service/TrackSystemSvc
    Utilities/DataHelper
    Utilities/KiTrack
)

set(SiliconTracking_srcs src/*.cpp)

# Modules
gaudi_add_module(SiliconTracking ${SiliconTracking_srcs}
    INCLUDE_DIRS GaudiKernel FWCore gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES TrackSystemSvcLib DataHelperLib KiTrackLib GaudiKernel FWCore ${GEAR_LIBRARIES} ${GSL_LIBRARIES} ${LCIO_LIBRARIES} 
)
