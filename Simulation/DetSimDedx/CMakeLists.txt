gaudi_subdir(DetSimDedx v0r0)

gaudi_depends_on_subdirs(
    FWCore
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimDedx_srcs
    src/DummyDedxSimTool.cpp
)

gaudi_add_module(DetSimDedx ${DetSimDedx_srcs}
    INCLUDE_DIRS
    LINK_LIBRARIES
        DD4hep
        ${DD4hep_COMPONENT_LIBRARIES}
        GaudiKernel
)

