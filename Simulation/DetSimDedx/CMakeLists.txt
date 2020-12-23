
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(EDM4HEP REQUIRED )
include_directories(${EDM4HEP_INCLUDE_DIR})

set(DetSimDedx_srcs
)

gaudi_add_module(DetSimDedx
                 SOURCES src/DummyDedxSimTool.cpp
                         src/BetheBlochEquationDedxSimTool.cpp
                         src/GFDndxSimTool.cpp

                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES}
                      Gaudi::GaudiKernel
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

