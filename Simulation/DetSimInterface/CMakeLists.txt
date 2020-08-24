
gaudi_subdir(DetSimInterface v0r0)

gaudi_depends_on_subdirs(GaudiKernel)

# DetSimInterface (headers only)
set(DetSimInterface_srcs
    src/IDetSimSvc.cpp
)

gaudi_add_library(DetSimInterface ${DetSimInterface_srcs}
    PUBLIC_HEADERS DetSimInterface
    LINK_LIBRARIES GaudiKernel
)

