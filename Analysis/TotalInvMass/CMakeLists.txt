gaudi_subdir(TotalInvMass v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(podio REQUIRED ) 
find_package(k4FWCore REQUIRED)

set(TotalInvMass_srcs src/TotalInvMass.cc)

gaudi_add_module(TotalInvMass ${TotalInvMass_srcs}
    INCLUDE_DIRS k4FWCore GaudiKernel GaudiAlgLib ${CLHEP_INCLUDE_DIR} DD4hep gear ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES k4FWCore GaudiKernel GaudiAlgLib ${CLHEP_LIBRARIES} DD4hep ${GEAR_LIBRARIES} ${GSL_LIBRARIES} ${LCIO_LIBRARIES}
    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
