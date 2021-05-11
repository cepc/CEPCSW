###############################################################################
# Package: DetCRD
#   CEPC Reference Detector (CRD)
################################################################################

find_package(DD4hep COMPONENTS DDRec DDG4 DDParsers REQUIRED)
# find_package(DD4hep)
find_package(Geant4)
include(${Geant4_USE_FILE})

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
include( DD4hep )

find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)

# install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION Detector/DetCRD)

gaudi_add_module(DetCRD
                 SOURCES src/Calorimeter/CRDEcal_v01.cpp
                         src/Other/CRDBeamPipe_v01_geo.cpp
                         src/Tracker/SiTrackerSkewRing_v01_geo.cpp

		 LINK ${DD4hep_COMPONENT_LIBRARIES}
)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(DetCRD)

install(TARGETS DetCRD
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
