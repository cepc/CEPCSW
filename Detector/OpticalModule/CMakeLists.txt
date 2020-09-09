###############################################################################
# Package: OpticalModule 
# modules for optical study  
################################################################################
gaudi_subdir(OpticalModule v0r0)

gaudi_depends_on_subdirs(GaudiKernel)

find_package(DD4hep COMPONENTS DDRec DDG4 DDParsers REQUIRED)
find_package(Geant4)
include(${Geant4_USE_FILE})

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
include( DD4hep )

find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/compact DESTINATION Detector/OpticalModule)

set(OpticalModule_src src/*.cpp)

gaudi_add_module(OpticalModule ${OpticalModule_src}
		 INCLUDE_DIRS
		 LINK_LIBRARIES 
                   DD4hep ${DD4hep_COMPONENT_LIBRARIES} 
)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(OpticalModule)
