gaudi_subdir(DCHDigi v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED )
include_directories(${EDM4HEP_INCLUDE_DIR})

find_package(podio REQUIRED )
find_package(ROOT COMPONENTS MathCore Physics GenVector Geom REQUIRED)

set(srcs
    src/*.cpp
)

gaudi_depends_on_subdirs(
    Detector/DetInterface
    Detector/DetSegmentation
)
## Modules
gaudi_add_module(DCHDigi ${srcs}
    INCLUDE_DIRS k4FWCore GaudiKernel GaudiAlgLib ${CLHEP_INCLUDE_DIR} DD4hep ROOT 
    LINK_LIBRARIES k4FWCore GaudiKernel GaudiAlgLib ${CLHEP_LIBRARIES} DD4hep ${DD4hep_COMPONENT_LIBRARIES} DDRec ROOT DetSegmentation
    -Wl,--no-as-needed 
    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
