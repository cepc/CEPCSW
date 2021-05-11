
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

gaudi_add_module(DetSimDedx
                 SOURCES src/DummyDedxSimTool.cpp
                         src/BetheBlochEquationDedxSimTool.cpp
                         src/GFDndxSimTool.cpp

                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES}
                      Gaudi::GaudiKernel
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

