gaudi_subdir(Tracking v0r0)

find_package(CLHEP REQUIRED)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(EDM4HEP REQUIRED ) 

gaudi_depends_on_subdirs(
    Service/GearSvc
    Service/EventSeeder
    Service/TrackSystemSvc
)

set(Tracking_srcs
    src/Clupatra/*.cpp
)

# Modules
gaudi_add_module(Tracking ${Tracking_srcs}
    INCLUDE_DIRS GaudiKernel FWCore CLHEP gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES GaudiKernel FWCore CLHEP $ENV{GEAR}/lib/libgearsurf.so ${GSL_LIBRARIES} ${LCIO_LIBRARIES} TrackSystemSvcLib
)
