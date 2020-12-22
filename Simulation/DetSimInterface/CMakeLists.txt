

# DetSimInterface (headers only)
set(DetSimInterface_srcs
    src/IDetSimSvc.cpp
)

gaudi_add_library(DetSimInterface
                 SOURCES ${DetSimInterface_srcs}
                 LINK Gaudi::GaudiKernel
)


