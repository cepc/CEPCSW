# Headers and Libraries

# Modules
gaudi_add_module(VXDTestbeam
                 SOURCES src/TBDataReader/TBDataReader.cpp
                 SOURCES src/TBDataReader/TBDataClusterReader.cpp
                 LINK #DetInterface
                      k4FWCore::k4FWCore
                      Gaudi::GaudiAlgLib Gaudi::GaudiKernel
                      ${LCIO_LIBRARIES} 
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ${podio_LIBRARIES} podio::podioRootIO
)

if (GenFit_FOUND)
  target_link_libraries(VXDTestbeam PUBLIC GenFit::genfit2)
endif()

target_include_directories(VXDTestbeam PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS VXDTestbeam
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
