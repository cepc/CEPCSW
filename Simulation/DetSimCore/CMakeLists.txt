
gaudi_subdir(DetSimCore v0r0)

gaudi_depends_on_subdirs(
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

set(DetSimCore_srcs
    src/DetSimAlg.cpp
    src/DetSimSvc.cpp
    src/DetectorConstruction.cpp
    src/PrimaryGeneratorAction.cpp
    src/ActionInitialization.cpp
    src/RunAction.cpp
    src/EventAction.cpp
)

gaudi_add_module(DetSimCore ${DetSimCore_srcs}
    INCLUDE_DIRS DetSimInterface GaudiKernel Geant4
    LINK_LIBRARIES DetSimInterface GaudiKernel Geant4
)
