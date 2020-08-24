##############################################################################
# Package: KalDet
#    Desc: import from ILCSoft
##############################################################################

gaudi_subdir(KalDet v0r0)

find_package(CLHEP REQUIRED)
find_package(LCIO)
find_package(GEAR)
find_package(ROOT COMPONENTS MathCore)
find_package(EDM4HEP REQUIRED)

gaudi_depends_on_subdirs(
    Detector/DetInterface
    Utilities/KalTest
)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

INCLUDE( MyFindROOT )
INCLUDE( MacroRootDict )

SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )

INCLUDE( ${ROOT_DICT_MACROS_FILE} )

SET( lib_input_dirs src/gen src/kern src/lctpc/gearTPC )

FOREACH( lib_input_dir ${lib_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${lib_input_dir} )
ENDFOREACH()

LIST( APPEND ROOT_DICT_INCLUDE_DIRS /workfs/bes/fucd/Key4hep/CEPCSW/Utilities/KalTest )

#MESSAGE( STATUS "ROOT_DICT_INCLUDE_DIRS: ${ROOT_DICT_INCLUDE_DIRS}" )

FOREACH( lib_input_dir ${lib_input_dirs} )
    AUX_SOURCE_DIRECTORY( ${lib_input_dir} lib_sources )
    PREPARE_ROOT_DICT_HEADERS( ${lib_input_dir} )
    INSTALL(DIRECTORY ${lib_input_dir}/ DESTINATION "include/kaldet"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE
    )
    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )
    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )
ENDFOREACH()

FILE( GLOB_RECURSE ild_sources "src/ild/*.cc" )
FILE( GLOB_RECURSE lib_headers "src/ild/*.h" )

include_directories(src/ild/common)

set(KalDetLib_srcs ${lib_sources} src/ild/*/*.cc )

gaudi_add_library(KalDetLib ${KalDetLib_srcs}
		 PUBLIC_HEADERS kaldet
                 LINK_LIBRARIES GaudiKernel ROOT CLHEP LCIO $ENV{GEAR}/lib/libgearsurf.so KalTestLib EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
