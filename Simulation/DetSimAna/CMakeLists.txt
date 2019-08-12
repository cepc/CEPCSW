
gaudi_subdir(DetSimAna v0r0)

gaudi_depends_on_subdirs(
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

set(DetSimAna_srcs
    src/ExampleAnaElemTool.cpp
)

gaudi_add_module(DetSimAna ${DetSimAna_srcs}
    INCLUDE_DIRS DetSimInterface GaudiKernel Geant4
    LINK_LIBRARIES DetSimInterface GaudiKernel Geant4
)
