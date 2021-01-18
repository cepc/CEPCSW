###############################################################################
# Package: DetCEPCv4 
# Ref to Package: DetFCCeeIDEA
# Based on package: lcgeo
################################################################################

find_package(DD4hep COMPONENTS DDRec DDG4 DDParsers REQUIRED)
# find_package(DD4hep)
find_package(Geant4)
include(${Geant4_USE_FILE})

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
include( DD4hep )

find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)

# install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION Detector/DetCEPCv4)

gaudi_add_module(DetCEPCv4
                 SOURCES src/tracker/VXD04_geo.cpp
                         src/tracker/FTD_Simple_Staggered_geo.cpp
                         src/tracker/FTD_cepc_geo.cpp
                         src/tracker/SIT_Simple_Pixel_geo.cpp 
                         src/tracker/SIT_Simple_Planar_geo.cpp
                         src/tracker/TPC10_geo.cpp
                         src/tracker/SET_Simple_Planar_geo.cpp
                         src/calorimeter/SEcal05_Helpers.cpp
                         src/calorimeter/SEcal05_Barrel.cpp
                         src/calorimeter/SEcal05_Endcaps.cpp
                         src/calorimeter/SEcal05_ECRing.cpp
                         src/calorimeter/SHcalRpc01_Barrel.cpp
                         src/calorimeter/SHcalRpc01_Endcaps.cpp
                         src/calorimeter/SHcalRpc01_EndcapRing.cpp
			 src/calorimeter/Yoke05_Barrel.cpp
			 src/calorimeter/Yoke05_Endcaps.cpp
                         src/other/BoxSupport_o1_v01_geo.cpp
                         src/other/TubeSupport_o1_v01_geo.cpp
			 src/other/SCoil02_geo.cpp

		 LINK ${DD4hep_COMPONENT_LIBRARIES}
)

target_include_directories(DetCEPCv4 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/src/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(DetCEPCv4)

install(TARGETS DetCEPCv4
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
