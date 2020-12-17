gaudi_subdir(DCHDedx v0r0)

find_package(EDM4HEP REQUIRED ) 
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)


gaudi_depends_on_subdirs(
    k4FWCore
    Simulation/DetSimInterface
    Detector/DetInterface
    Detector/DetSegmentation
)

set(srcs
    src/*.cpp
)

# Modules
gaudi_add_module(DCHDedx ${srcs}
    INCLUDE_DIRS GaudiKernel 
    LINK_LIBRARIES GaudiAlgLib GaudiKernel 
    DetSegmentation
      -Wl,--no-as-needed
     EDM4HEP::edm4hep EDM4HEP::edm4hepDict
     -Wl,--as-needed
)
