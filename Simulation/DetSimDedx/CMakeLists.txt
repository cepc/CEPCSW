gaudi_subdir(DetSimDedx v0r0)

gaudi_depends_on_subdirs(
    k4FWCore
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED )
include_directories(${EDM4HEP_INCLUDE_DIR})

set(DetSimDedx_srcs
    src/DummyDedxSimTool.cpp
    src/BetheBlochEquationDedxSimTool.cpp
)

gaudi_add_module(DetSimDedx ${DetSimDedx_srcs}
    INCLUDE_DIRS
    LINK_LIBRARIES
        DD4hep
        ${DD4hep_COMPONENT_LIBRARIES}
        GaudiKernel
    -Wl,--no-as-needed 
    EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

