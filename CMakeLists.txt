CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

find_package(ROOT COMPONENTS RIO Tree)

find_package(GaudiProject)

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "14|17")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()


gaudi_project(CEPCSW v0r1
              USE Gaudi v33r1 
              USE K4FWCore v0r1
)