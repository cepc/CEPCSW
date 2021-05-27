#################################################################################
##Package : MagneticFieldMap
#################################################################################

gaudi_add_module(MagneticFieldMap
                 SOURCES src/GenericBFieldMapBrBz.cpp
                         src/GenericBFieldMapBrBzFactory.cpp
                         src/FieldMapFileProvider.cpp
                 LINK Gaudi::GaudiKernel
                      ${DD4hep_COMPONENT_LIBRARIES}
                      ${ROOT_LIBRARIES}
)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(MagneticFieldMap)

install(TARGETS MagneticFieldMap
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
