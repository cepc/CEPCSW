gaudi_subdir(Tracking v0r0)

find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(EDM4HEP REQUIRED ) 
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)


gaudi_depends_on_subdirs(
    Service/GearSvc
    Service/EventSeeder
    Service/TrackSystemSvc
    Detector/DetSegmentation
)

set(Tracking_srcs
    src/Clupatra/*.cpp
    src/FullLDCTracking/*.cpp
    src/TruthTracker/*.cpp
)

# Modules
gaudi_add_module(Tracking ${Tracking_srcs}
    INCLUDE_DIRS GaudiKernel gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES GaudiAlgLib GaudiKernel ${GEAR_LIBRARIES} ${GSL_LIBRARIES} ${LCIO_LIBRARIES} TrackSystemSvcLib
    DetSegmentation
      -Wl,--no-as-needed
     EDM4HEP::edm4hep EDM4HEP::edm4hepDict
     -Wl,--as-needed
)
