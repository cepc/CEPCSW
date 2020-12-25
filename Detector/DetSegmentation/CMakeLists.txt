#################################################################################
##Package : DetSegmentation
#################################################################################

find_package(DD4hep COMPONENTS DDRec DDG4 DDParsers REQUIRED)
set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )

include(DD4hep)

gaudi_add_library(DetSegmentation
                  SOURCES src/GridDriftChamber.cpp
                  LINK Gaudi::GaudiKernel
                       ${DD4hep_COMPONENT_LIBRARIES}
                       ${ROOT_LIBRARIES}
)

gaudi_add_module(DetSegmentationPlugin
                 SOURCES src/plugins/SegmentationFactories.cpp
                 LINK Gaudi::GaudiKernel
                      ${ROOT_LIBRARIES}
                      ${DD4hep_COMPONENT_LIBRARIES}
                      DetSegmentation)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(DetSegmentationPlugin)

install(TARGETS DetSegmentation DetSegmentationPlugin
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
