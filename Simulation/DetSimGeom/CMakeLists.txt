
gaudi_subdir(DetSimGeom v0r0)

gaudi_depends_on_subdirs(
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

set(DetSimGeom_srcs
    src/WorldDetElemTool.cpp
)

gaudi_add_module(DetSimGeom ${DetSimGeom_srcs}
    INCLUDE_DIRS DetSimInterface GaudiKernel Geant4
    LINK_LIBRARIES DetSimInterface GaudiKernel Geant4
)
