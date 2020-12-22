
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimFastModel_srcs
    src/DummyFastSimG4Tool.cpp
    src/DummyFastSimG4Model.cpp
)

gaudi_add_module(DetSimFastModel
                 SOURCES ${DetSimFastModel_srcs}
                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES}
                      Gaudi::GaudiKernel
)

