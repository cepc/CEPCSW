
gaudi_subdir(Examples v0r0)

set(Examples_srcs
    src/HelloWorld/*.cpp
    src/FirstSvc/*.cpp
    src/SecondAlg/*.cpp
)

# Headers and Libraries

gaudi_install_headers(Examples)

# Modules
gaudi_add_module(Examples ${Examples_srcs}
    INCLUDE_DIRS GaudiKernel
    LINK_LIBRARIES GaudiKernel
)

# Unit tests
gaudi_add_test(HelloAlg
               FRAMEWORK options/helloalg.py)