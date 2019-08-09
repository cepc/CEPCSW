
gaudi_subdir(DetSimInterface v0r0)

# DetSimInterface (headers only)
set(DetSimInterface_srcs
    src/*.cpp
)

gaudi_add_library(DetSimInterface ${DetSimInterface_srcs}
    PUBLIC_HEADERS DetSimInterface
    LINK_LIBRARIES GaudiKernel
)

