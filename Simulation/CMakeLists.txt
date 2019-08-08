
gaudi_subdir(Simulation v0r0)

find_package(Geant4)
include(${Geant4_USE_FILE})

# DetSimCore module

set(DetSimCore_srcs
    src/DetSimCore/*.cpp
)

gaudi_add_module(DetSimCore ${DetSimCore_srcs}
    INCLUDE_DIRS GaudiKernel Geant4
    LINK_LIBRARIES GaudiKernel Geant4
)
