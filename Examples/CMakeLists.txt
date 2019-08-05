
gaudi_subdir(Examples v0r0)

set(Examples_srcs
    src/HelloWorld/*.cpp
)

gaudi_add_module(Examples ${Examples_srcs}
    INCLUDE_DIRS GaudiKernel
    LINK_LIBRARIES GaudiKernel
)