##############################################################################
# Package: KalTest
#    Desc: import from ILCSoft
##############################################################################

gaudi_subdir(KalTest v0r0)

find_package(ROOT REQUIRED COMPONENTS MathCore)

gaudi_depends_on_subdirs()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )

INCLUDE( MyFindROOT )
INCLUDE( MacroRootDict )

SET( lib_input_dirs src/geomlib src/kallib src/kaltracklib src/utils )

FOREACH( lib_input_dir ${lib_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${lib_input_dir} )
ENDFOREACH()

#MESSAGE( STATUS "ROOT_DICT_INCLUDE_DIRS: ${ROOT_DICT_INCLUDE_DIRS}" )

FOREACH( lib_input_dir ${lib_input_dirs} )

    AUX_SOURCE_DIRECTORY( ${lib_input_dir} lib_sources )

    PREPARE_ROOT_DICT_HEADERS( ${lib_input_dir} )

    INSTALL(DIRECTORY ${lib_input_dir}/ DESTINATION "include/kaltest"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE
    )

    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )

    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

include_directories( ${ROOT_DICT_INCLUDE_DIRS} )
include_directories( ${ROOT_INCLUDE_DIRS} )

set(KalTestLib_srcs ${lib_sources})

gaudi_add_library(KalTestLib ${KalTestLib_srcs}
			     PUBLIC_HEADERS kaltest
			     LINK_LIBRARIES GaudiKernel ROOT
)
