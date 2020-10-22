gaudi_subdir(G2CDArbor v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(podio REQUIRED ) 
find_package(K4FWCore REQUIRED)

message("EDM4HEP_INCLUDE_DIRS: ${EDM4HEP_INCLUDE_DIR}")
message("EDM4HEP_LIB: ${EDM4HEP_LIBRARIES}")
include_directories(${EDM4HEP_INCLUDE_DIR})

gaudi_depends_on_subdirs(
    Service/GearSvc
    Detector/DetInterface
)

set(G2CDArbor_srcs src/*.cpp)

# Modules
gaudi_add_module(G2CDArbor ${G2CDArbor_srcs}
    INCLUDE_DIRS K4FWCore GaudiKernel GaudiAlgLib ${CLHEP_INCLUDE_DIR} DD4hep gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES K4FWCore GaudiKernel GaudiAlgLib ${CLHEP_LIBRARIES} DD4hep ${GEAR_LIBRARIES} ${GSL_LIBRARIES} ${LCIO_LIBRARIES}
    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
