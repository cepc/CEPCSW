
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

set(DetSimSD_srcs
    src/CalorimeterSensDetTool.cpp

    src/DDG4SensitiveDetector.cpp
    src/CaloSensitiveDetector.cpp

    src/DriftChamberSensDetTool.cpp
    src/DriftChamberSensitiveDetector.cpp

    src/TimeProjectionChamberSensDetTool.cpp
    src/TimeProjectionChamberSensitiveDetector.cpp
)

gaudi_add_module(DetSimSD 
                 SOURCES ${DetSimSD_srcs}
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

