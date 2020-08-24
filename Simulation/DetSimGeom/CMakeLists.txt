
gaudi_subdir(DetSimGeom v0r0)

gaudi_depends_on_subdirs(
    Simulation/DetSimInterface
    Detector/DetInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimGeom_srcs
    src/WorldDetElemTool.cpp
    src/AnExampleDetElemTool.cpp
)

gaudi_add_module(DetSimGeom ${DetSimGeom_srcs}
    INCLUDE_DIRS
    # DetSimInterface
    # DetInterface
    # DD4hep
    # GaudiKernel
    # Geant4
    LINK_LIBRARIES
      DD4hep ${DD4hep_COMPONENT_LIBRARIES} 
      GaudiKernel 
      # Geant4
)
