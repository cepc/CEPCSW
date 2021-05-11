################################################################################
# Package: Identifier
################################################################################

gaudi_add_header_only_library(Identifier)

install(TARGETS Identifier
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
