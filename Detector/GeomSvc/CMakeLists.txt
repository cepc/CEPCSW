
##############################################################################
# Package: GeomSvc
#    Desc: implement the IGeomSvc interface.
##############################################################################

find_package(ROOT COMPONENTS MathCore GenVector Geom Tree)
find_package(DD4hep COMPONENTS DDG4 DDRec REQUIRED)
find_package(GEAR REQUIRED)

gaudi_add_module(GeomSvc
                 SOURCES src/GeomSvc.cpp
                 LINK
                   DetInterface
                   ${DD4hep_COMPONENT_LIBRARIES} 
                   Gaudi::GaudiKernel
                   ${GEAR_LIBRARIES}
                   ${ROOT_LIBRARIES}
)

install(TARGETS GeomSvc
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
