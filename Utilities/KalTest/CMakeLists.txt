##############################################################################
# Package: KalTest
#    Desc: import from ILCSoft
##############################################################################

gaudi_subdir(KalTest v0r0)

find_package(ROOT REQUIRED COMPONENTS MathCore)

gaudi_depends_on_subdirs()

set( DICT_CINT_DEFINITIONS "HANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
set( DICT_INPUT_DIRS geomlib kallib kaltracklib utils )
set( DICT_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rootdict )

foreach( DICT ${DICT_INPUT_DIRS} )
  list( APPEND DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT} )
  aux_source_directory( src/${DICT} LIB_SOURCES )
endforeach()

foreach( DICT ${DICT_INPUT_DIRS} )
  file( GLOB DICT_HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT}/*.h" )
  file( GLOB LINKDEFS "${CMAKE_CURRENT_SOURCE_DIR}/src/${DICT}/LinkDef.h" )

  if( LINKDEFS )
    list( REMOVE_ITEM DICT_HEADER_FILES "${LINKDEFS}" )
  endif()
    
  add_custom_command(OUTPUT ${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm
                     COMMAND ${ROOT_rootcling_CMD}
                     ARGS -f ${DICT_OUTPUT_DIR}/${DICT}.cxx -s ${DICT_OUTPUT_DIR}/${DICT} -inlineInputHeader
                     "$<$<BOOL:$<JOIN:${DICT_CINT_DEFINITIONS},>>:-D$<JOIN:${DICT_CINT_DEFINITIONS},;-D>>"
                     "$<$<BOOL:$<JOIN:${DICT_INCLUDE_DIRS},>>:-I$<JOIN:${DICT_INCLUDE_DIRS},;-I>>"
                     "$<JOIN:${DICT_HEADER_FILES},;>" "$<JOIN:${LINKDEFS},;>"
                     DEPENDS ${DICT_HEADER_FILES} ${LINKDEFS}
                     COMMAND_EXPAND_LISTS
                     )
  add_custom_target( ${DICT} DEPENDS ${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm )

  set_source_files_properties(${DICT_OUTPUT_DIR}/${DICT}.cxx ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm
                              PROPERTIES
                              GENERATED TRUE
                              COMPILE_FLAGS "-Wno-unused-function -Wno-overlength-strings"
			      )

  install(FILES ${DICT_OUTPUT_DIR}/${DICT}_rdict.pcm DESTINATION lib)

  list( APPEND LIB_SOURCES ${DICT_OUTPUT_DIR}/${DICT}.cxx )
endforeach()

include_directories( ${DICT_INCLUDE_DIRS} )
include_directories( ${ROOT_INCLUDE_DIRS} )

set(KalTestLib_srcs ${LIB_SOURCES})

gaudi_add_library(KalTestLib ${KalTestLib_srcs}
			     PUBLIC_HEADERS kaltest
			     LINK_LIBRARIES GaudiKernel ROOT
)
