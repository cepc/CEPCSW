
gaudi_subdir(DetSimAna v0r0)

gaudi_depends_on_subdirs(
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimAna_srcs
    src/ExampleAnaElemTool.cpp
)

gaudi_add_module(DetSimAna ${DetSimAna_srcs}
    INCLUDE_DIRS DetSimInterface DD4hep GaudiKernel Geant4
    LINK_LIBRARIES DetSimInterface DD4hep ${DD4hep_COMPONENT_LIBRARIES} GaudiKernel Geant4
)
