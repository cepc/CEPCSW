gaudi_subdir(SimpleDigi v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(podio REQUIRED ) 
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)

gaudi_depends_on_subdirs(
    Service/GearSvc
    Service/EventSeeder
    Service/TrackSystemSvc
    Utilities/DataHelper
)

set(SimpleDigi_srcs src/*.cpp)

# Modules
gaudi_add_module(SimpleDigi ${SimpleDigi_srcs}
    INCLUDE_DIRS k4FWCore GaudiKernel GaudiAlgLib ${CLHEP_INCLUDE_DIR} gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES k4FWCore GaudiKernel GaudiAlgLib ${CLHEP_LIBRARIES} ${GEAR_LIBRARIES} ${GSL_LIBRARIES} ${LCIO_LIBRARIES} EDM4HEP::edm4hep EDM4HEP::edm4hepDict DataHelperLib
)
