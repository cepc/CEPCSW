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

if(ROOT_FOUND)
    message("found ROOT: ${ROOT_INCLUDE_DIRS} ${ROOT_LIBRARIES}")
endif(ROOT_FOUND)
if(LCIO_FOUND)
    message("found LCIO: ${LCIO_INCLUDE_DIRS} ${LCIO_LIBRARIES}")
endif(LCIO_FOUND)
if(podio_FOUND)
    message("found podio: ${podio_INCLUDE_DIRS} ${podio_LIBRARIES}")
endif(podio_FOUND)
if(HepMC_FOUND)
    message("found HepMC: ${HepMC_INCLUDE_DIRS} ${HepMC_LIBRARY_DIR}")
endif(HepMC_FOUND)
if(CLHEP_FOUND)
    message("found CLHEP: ${CLHEP_INCLUDE_DIRS} ${CLHEP_LIBRARY_DIR}")
endif(CLHEP_FOUND)

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
  )

install(TARGETS GenAlgo
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
