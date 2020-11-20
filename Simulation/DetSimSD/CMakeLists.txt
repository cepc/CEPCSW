
gaudi_subdir(DetSimSD v0r0)

gaudi_depends_on_subdirs(
    k4FWCore
    Detector/DetInterface
    Simulation/DetSimInterface
)

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

gaudi_add_module(DetSimSD ${DetSimSD_srcs}
    INCLUDE_DIRS
      # DetSimInterface
      # k4FWCore
      # DD4hep
      # GaudiKernel
      # Geant4 
    LINK_LIBRARIES
      # DetSimInterface
      # k4FWCore
      DD4hep
      ${DD4hep_COMPONENT_LIBRARIES} 
      GaudiKernel
      # Geant4
)
