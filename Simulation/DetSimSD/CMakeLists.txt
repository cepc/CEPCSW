
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

gaudi_add_module(DetSimSD 
                 SOURCES src/CalorimeterSensDetTool.cpp
                         src/DDG4SensitiveDetector.cpp
                         src/CaloSensitiveDetector.cpp

                         src/DriftChamberSensDetTool.cpp
                         src/DriftChamberSensitiveDetector.cpp

                         src/TimeProjectionChamberSensDetTool.cpp
                         src/TimeProjectionChamberSensitiveDetector.cpp

                 LINK DetSimInterface
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

