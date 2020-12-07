gaudi_subdir(DetSimFastModel v0r0)

gaudi_depends_on_subdirs(
    k4FWCore
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimFastModel_srcs
    src/DummyFastSimG4Tool.cpp
    src/DummyFastSimG4Model.cpp
)

gaudi_add_module(DetSimFastModel ${DetSimFastModel_srcs}
    INCLUDE_DIRS
    LINK_LIBRARIES
        DD4hep
        ${DD4hep_COMPONENT_LIBRARIES}
        GaudiKernel
)

