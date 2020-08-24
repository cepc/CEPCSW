###############################################################################
# Package: DetCEPCv4 
# Ref to Package: DetFCCeeIDEA
# Based on package: lcgeo
################################################################################
gaudi_subdir(DetCEPCv4 v0r0)

gaudi_depends_on_subdirs(GaudiKernel)


find_package(DD4hep COMPONENTS DDRec DDG4 DDParsers REQUIRED)
# find_package(DD4hep)
find_package(Geant4)
include(${Geant4_USE_FILE})

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
include( DD4hep )

find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION Detector/DetCEPCv4)

set(DetCEPCv4_src
    src/tracker/VXD04_geo.cpp
    src/tracker/FTD_Simple_Staggered_geo.cpp
    src/tracker/SIT_Simple_Pixel_geo.cpp 
    src/tracker/TPC10_geo.cpp
    src/tracker/SET_Simple_Planar_geo.cpp
    src/calorimeter/SEcal05_Helpers.cpp
    src/calorimeter/SEcal05_Barrel.cpp
    src/calorimeter/SEcal05_Endcaps.cpp
    src/calorimeter/SEcal05_ECRing.cpp
    src/other/BoxSupport_o1_v01_geo.cpp
    src/other/TubeSupport_o1_v01_geo.cpp
)

gaudi_add_module(DetCEPCv4
                  ${DetCEPCv4_src}
		  INCLUDE_DIRS 
                  # DD4hep 
                  # ROOT 
                  # Geant4 
                  src/include
		  LINK_LIBRARIES
                  # GaudiKernel
                  DD4hep
                  ${DD4hep_COMPONENT_LIBRARIES}
                  # ROOT 
                  # Geant4
)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(DetCEPCv4)
