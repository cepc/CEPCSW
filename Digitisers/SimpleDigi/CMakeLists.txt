gaudi_subdir(SimpleDigi v0r0)

find_package(CLHEP REQUIRED)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(podio REQUIRED ) 
find_package(K4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)

gaudi_depends_on_subdirs(
    Service/GearSvc
    Service/EventSeeder
    Service/TrackSystemSvc
)

set(SimpleDigi_srcs src/*.cpp)

# Modules
gaudi_add_module(SimpleDigi ${SimpleDigi_srcs}
    INCLUDE_DIRS K4FWCore GaudiKernel GaudiAlgLib CLHEP gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES K4FWCore GaudiKernel GaudiAlgLib CLHEP $ENV{GEAR}/lib/libgearsurf.so ${GSL_LIBRARIES} ${LCIO_LIBRARIES} EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
