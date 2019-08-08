
gaudi_subdir(Simulation v0r0)

find_package(Geant4)
include(${Geant4_USE_FILE})

# DetSimInterface (headers only)
set(DetSimInterface_srcs
    src/DetSimInterface/*.cpp
)


gaudi_add_library(DetSimInterfaces ${DetSimInterface_srcs}
    INCLUDE_DIRS include
    PUBLIC_HEADERS include/DetSimInterface
)

# DetSimCore module

set(DetSimCore_srcs
    src/DetSimCore/*.cpp
)

gaudi_add_module(DetSimCore ${DetSimCore_srcs}
    INCLUDE_DIRS DetSimInterfaces GaudiKernel Geant4
    LINK_LIBRARIES DetSimInterfaces GaudiKernel Geant4
)
