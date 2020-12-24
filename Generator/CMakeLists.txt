########################################

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT COMPONENTS RIO Tree TreePlayer MathCore Net Graf3d Graf Gpad EG REQUIRED)
find_package(LCIO)
find_package(podio)
find_package(EDM4HEP)
find_package(HepMC)
find_package(CLHEP REQUIRED;CONFIG)
find_package(k4FWCore REQUIRED)

gaudi_add_module(GenAlgo
                 SOURCES src/IGenTool.cpp 
                         src/GenAlgo.cpp 
                         src/GenEvent.cpp 
                         src/GenReader.cpp 
                         src/StdHepRdr.cpp 
                         src/GenPrinter.cpp
                         # src/LCAscHepRdr.cc
                         # src/HepevtRdr.cpp
                         src/SLCIORdr.cpp
                         src/HepMCRdr.cpp
                         src/GtGunTool.cpp

                 LINK ${ROOT_LIBRARIES}
                      k4FWCore::k4FWCore 
                      Gaudi::GaudiAlgLib
                      Gaudi::GaudiKernel
                      ${HEPMC_LIBRARIES}
                      ${CLHEP_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ROOT::EG
)

install(TARGETS GenAlgo
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
