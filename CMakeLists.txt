CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

project(CEPCSW)

find_package(ROOT COMPONENTS RIO Tree)

find_package(Gaudi)

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "14|17")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

add_subdirectory(Analysis)
add_subdirectory(Detector)
