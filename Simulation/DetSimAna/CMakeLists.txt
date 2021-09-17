
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

gaudi_add_module(DetSimAna
                 SOURCES src/Edm4hepWriterAnaElemTool.cpp
                 LINK DetSimInterface
                      ${DD4hep_COMPONENT_LIBRARIES} 
                      Gaudi::GaudiKernel
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      k4FWCore::k4FWCore
)

install(TARGETS DetSimAna
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
