
gaudi_subdir(DetSimAna v0r0)

gaudi_depends_on_subdirs(
    FWCore
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

# For EDM & I/O
find_package(podio REQUIRED)
find_package(plcio REQUIRED)
# find_package(LCIO REQUIRED)

set(DetSimAna_srcs
    src/ExampleAnaElemTool.cpp
)
message("podio_LIBRARIES: ${podio_LIBRARIES}")
message("plcio_LIBRARIES: ${plcio_LIBRARIES}")
gaudi_add_module(DetSimAna ${DetSimAna_srcs}
    INCLUDE_DIRS DetSimInterface FWCore DD4hep GaudiKernel Geant4 
                 ${plcio_INCLUDE_DIRS} ${podio_INCLUDE_DIRS}
    LINK_LIBRARIES DetSimInterface FWCore DD4hep ${DD4hep_COMPONENT_LIBRARIES} GaudiKernel Geant4
                 ${plcio_LIBRARIES} ${podio_LIBRARIES}
)
