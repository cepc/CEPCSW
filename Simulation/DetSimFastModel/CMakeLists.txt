
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

gaudi_add_module(DetSimFastModel
                 SOURCES src/DummyFastSimG4Tool.cpp
                         src/DummyFastSimG4Model.cpp
                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES}
                      Gaudi::GaudiKernel
)

