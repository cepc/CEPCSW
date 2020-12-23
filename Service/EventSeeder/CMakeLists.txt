
set(EventSeeder_srcs
    src/MarlinEvtSeeder.cpp
)

gaudi_add_header_only_library(EventSeeder)

gaudi_add_module(EventSeederPlugins
                 SOURCES ${EventSeeder_srcs}
                 LINK EventSeeder
                      Gaudi::GaudiKernel
)
