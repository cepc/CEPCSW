
gaudi_add_header_only_library(GearSvc)

gaudi_add_module(GearSvcPlugins
                 SOURCES src/GearSvc.cpp
                 LINK GearSvc
                      Gaudi::GaudiKernel
                      ${GEAR_LIBRARIES}
                      ${DD4hep_COMPONENT_LIBRARIES}
                      DetInterface
		      DetSegmentation
)

install(TARGETS GearSvc GearSvcPlugins
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

