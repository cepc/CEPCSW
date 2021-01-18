CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

project(CEPCSW)

# setup some necessary envvar
include(cmake/CEPCSWEnv.cmake)
#

find_package(ROOT COMPONENTS RIO Tree)

find_package(Gaudi)


include(GNUInstallDirs)
include(CTest)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/InstallArea/ CACHE PATH
    "Install path prefix, prepended onto install directories." FORCE )
endif()

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "14|17")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake") # (Find*.cmake)
include(cmake/CEPCSWDependencies.cmake)

add_subdirectory(Analysis)
add_subdirectory(Detector)
add_subdirectory(Digitisers)
add_subdirectory(Examples)
add_subdirectory(Generator)
add_subdirectory(Reconstruction)
add_subdirectory(Service)
add_subdirectory(Simulation)
add_subdirectory(Utilities)

##############################################################################
# INSTALL
##############################################################################

install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/"
  )

gaudi_install(CMAKE cmake/${PROJECT_NAME}Config.cmake)
