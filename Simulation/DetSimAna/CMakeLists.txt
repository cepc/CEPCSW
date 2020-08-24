
gaudi_subdir(DetSimAna v0r0)

gaudi_depends_on_subdirs(
    FWCore
    Simulation/DetSimInterface
)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

# For EDM & I/O
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
# find_package(LCIO REQUIRED)

set(DetSimAna_srcs
#    src/ExampleAnaElemTool.cpp
    src/Edm4hepWriterAnaElemTool.cpp
)
message("podio_LIBRARIES: ${podio_LIBRARIES}")
gaudi_add_module(DetSimAna ${DetSimAna_srcs}
    INCLUDE_DIRS
      # DetSimInterface
      # FWCore
      # DD4hep
      # GaudiKernel
      # Geant4 
      # ${plcio_INCLUDE_DIRS}
      # ${podio_INCLUDE_DIRS}
    LINK_LIBRARIES
      # DetSimInterface
      # FWCore
      DD4hep
      ${DD4hep_COMPONENT_LIBRARIES} 
      GaudiKernel
      # Geant4

      -Wl,--no-as-needed 
      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
      -Wl,--as-needed 
      # ${podio_LIBRARIES}
)
