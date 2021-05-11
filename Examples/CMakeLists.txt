# Headers and Libraries

# Modules
gaudi_add_module(Examples
                 SOURCES src/HelloWorld/HelloAlg.cpp
                         src/FirstSvc/FirstSvc.cpp
                         src/SecondAlg/SecondAlg.cpp
                         src/Edm4hepTest/Edm4hepReadAlg.cpp
                         src/Edm4hepTest/Edm4hepReadDCAlg.cpp
                         src/Edm4hepTest/Edm4hepWriteAlg.cpp
                         src/DumpIDAlg/DumpIDAlg.cpp
                 LINK DetInterface
                      k4FWCore::k4FWCore
                      Gaudi::GaudiAlgLib Gaudi::GaudiKernel
                      ${LCIO_LIBRARIES} 
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ${podio_LIBRARIES} podio::podioRootIO
)

if (GenFit_FOUND)
  target_link_libraries(Examples PUBLIC GenFit::genfit2)
endif()

target_include_directories(Examples PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS Examples
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
