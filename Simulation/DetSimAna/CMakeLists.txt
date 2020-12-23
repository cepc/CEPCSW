
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(DD4hep COMPONENTS DDG4 REQUIRED)

# For EDM & I/O
find_package(EDM4HEP REQUIRED)


gaudi_add_module(DetSimAna
                 SOURCES src/Edm4hepWriterAnaElemTool.cpp
                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES} 
                      Gaudi::GaudiKernel
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

install(TARGETS DetSimAna
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
