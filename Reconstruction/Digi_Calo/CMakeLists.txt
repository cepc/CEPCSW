gaudi_subdir(Digi_Calo v0r0)

find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED )
message("EDM4HEP_INCLUDE_DIRS: ${EDM4HEP_INCLUDE_DIR}")
message("EDM4HEP_LIB: ${EDM4HEP_LIBRARIES}")
include_directories(${EDM4HEP_INCLUDE_DIR})

find_package(CLHEP REQUIRED)
find_package(podio REQUIRED )

set(srcs
    src/*.cpp
)

gaudi_depends_on_subdirs(
    Detector/DetInterface
)
## Modules
gaudi_add_module(Digi_Calo ${srcs}
    INCLUDE_DIRS FWCore GaudiKernel GaudiAlgLib CLHEP DD4hep 
    LINK_LIBRARIES FWCore GaudiKernel GaudiAlgLib CLHEP DD4hep ${DD4hep_COMPONENT_LIBRARIES} DDRec
    -Wl,--no-as-needed 
    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
