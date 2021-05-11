# - Try to find GenFit
# Defines:
#
#  GenFit_FOUND
#  GenFit_INCLUDE_DIRS
#  GenFit_LIBRARIES
#
# Imports:
#
#  GenFit::genfit2
#
# Usage of the target instead of the variables is advised

# Find quietly if already found before
if(DEFINED CACHE{GenFit_INCLUDE_DIRS})
  set(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY YES)
endif()

find_path(GenFit_INCLUDE_DIRS NAMES GFGbl.h
          HINTS ${GENFIT_ROOT_DIR}/include $ENV{GENFIT_ROOT_DIR}/include
                ${GenFit_DIR}/include $ENV{GenFit_DIR}/include
                ${CMAKE_PREFIX_PATH}/include
)
find_library(GenFit_LIBRARIES NAMES genfit2
          HINTS $ENV{GENFIT_ROOT_DIR}/lib ${GENFIT_ROOT_DIR}/lib
                ${GenFit_DIR}/lib ${GenFit_DIR}/lib64
                $ENV{GenFit_DIR}/lib $ENV{GenFit_DIR}/lib64
                ${CMAKE_PREFIX_PATH}/lib ${CMAKE_PREFIX_PATH}/lib64
)

# handle the QUIETLY and REQUIRED arguments and set GENFIT_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GenFit DEFAULT_MSG GenFit_INCLUDE_DIRS GenFit_LIBRARIES)

mark_as_advanced(GenFit_FOUND GenFit_INCLUDE_DIRS GenFit_LIBRARIES)

# Modernisation: create an interface target to link against
if(TARGET GenFit::genfit2)
    return()
endif()
if(GenFit_FOUND)
  add_library(GenFit::genfit2 IMPORTED INTERFACE)
  target_include_directories(GenFit::genfit2 SYSTEM INTERFACE "${GenFit_INCLUDE_DIRS}")
  target_link_libraries(GenFit::genfit2 INTERFACE "${GenFit_LIBRARIES}")

  # Display the imported target for the user to know
  if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message(STATUS "  Import target: GenFit::genfit2")
  endif()
endif()

if(COMMAND __deprecate_var_for_target)
  variable_watch(GenFit_INCLUDE_DIRS __deprecate_var_for_target)
endif()
