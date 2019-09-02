gaudi_subdir(EventSeeder v0r0)

#find_package(GEAR REQUIRED)

set(EventSeeder_srcs
    src/*.cpp
)

gaudi_install_headers(EventSeeder)

gaudi_add_module(EventSeeder ${EventSeeder_srcs}
    INCLUDE_DIRS GaudiKernel
    LINK_LIBRARIES GaudiKernel
)
