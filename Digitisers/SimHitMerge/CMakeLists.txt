gaudi_subdir(SimHitMergeAlg v0r0)

find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(podio REQUIRED ) 
find_package(K4FWCore REQUIRED)

include_directories(${EDM4HEP_INCLUDE_DIR})

gaudi_depends_on_subdirs(
    Detector/DetInterface
)
set(SimHitMergeAlg_srcs src/*.cpp)

# Modules
gaudi_add_module(SimHitMerge ${SimHitMergeAlg_srcs}
    INCLUDE_DIRS K4FWCore GaudiKernel GaudiAlgLib DD4hep  
    LINK_LIBRARIES K4FWCore GaudiKernel GaudiAlgLib DD4hep DDRec
    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
