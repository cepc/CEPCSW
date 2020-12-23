
gaudi_add_header_only_library(EventSeeder)

gaudi_add_module(EventSeederPlugins
                 SOURCES src/MarlinEvtSeeder.cpp
                 LINK EventSeeder
                      Gaudi::GaudiKernel
)

install(TARGETS EventSeeder EventSeederPlugins
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


