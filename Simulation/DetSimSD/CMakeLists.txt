
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

gaudi_add_library(DetSimSDLib
                  SOURCES src/Geant4Hits.cpp
                          src/DDG4SensitiveDetector.cpp
                          src/CaloSensitiveDetector.cpp
                          src/DriftChamberSensitiveDetector.cpp
                          src/TimeProjectionChamberSensitiveDetector.cpp
                          src/GenericTrackerSensitiveDetector.cpp
                          src/TrackerCombineSensitiveDetector.cpp
                  LINK DetSimInterface
                       DetInterface
                       ${DD4hep_COMPONENT_LIBRARIES}
)
target_include_directories(DetSimSDLib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS DetSimSDLib
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

gaudi_add_module(DetSimSD 
                 SOURCES src/CalorimeterSensDetTool.cpp
                         src/DriftChamberSensDetTool.cpp
                         src/TimeProjectionChamberSensDetTool.cpp
                         src/GenericTrackerSensDetTool.cpp
                 LINK DetSimInterface
                      DetSimSDLib
                      DetInterface
                      ${DD4hep_COMPONENT_LIBRARIES} 
                      Gaudi::GaudiKernel
)
target_include_directories(DetSimSD PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS DetSimSD
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

