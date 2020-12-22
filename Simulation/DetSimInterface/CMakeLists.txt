

# DetSimInterface (headers only)
set(DetSimInterface_srcs
    src/IDetSimSvc.cpp
)

gaudi_add_library(DetSimInterface
                 SOURCES ${DetSimInterface_srcs}
                 LINK Gaudi::GaudiKernel
)

install(TARGETS DetSimInterface
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

