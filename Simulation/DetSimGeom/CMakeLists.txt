
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimGeom_srcs
    src/WorldDetElemTool.cpp
    src/AnExampleDetElemTool.cpp
)

gaudi_add_module(DetSimGeom 
                 SOURCES ${DetSimGeom_srcs}
                 LINK DetSimInterface
                      DetInterface
                      ${DD4hep_COMPONENT_LIBRARIES} 
                      Gaudi::GaudiKernel 
                      ${Geant4_LIBRARIES}
)

install(TARGETS DetSimGeom
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)


