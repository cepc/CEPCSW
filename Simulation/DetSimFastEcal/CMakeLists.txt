
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})


gaudi_add_module(DetSimFastEcal
                 SOURCES src/EcalFastSimG4Tool.cpp
                         src/EcalFastSimG4Model.cpp
                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES}
                      Gaudi::GaudiKernel
)

