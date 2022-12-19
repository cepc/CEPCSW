

# DetSimInterface (headers only)

gaudi_add_library(DetSimInterface
                 SOURCES src/IDetSimSvc.cpp
                         src/CommonUserEventInfo.cc
                 LINK Gaudi::GaudiKernel
                      ${Geant4_LIBRARIES}
)

install(TARGETS DetSimInterface
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

