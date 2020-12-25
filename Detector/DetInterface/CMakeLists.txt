################################################################################
# Package: DetInterface
################################################################################

gaudi_add_header_only_library(DetInterface)

install(TARGETS DetInterface
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
