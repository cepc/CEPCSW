########################################
gaudi_subdir(Generator v0r0)
set(GenAlgo_srcs 
    src/IGenTool.cpp 
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
)
set(GenAlgo_incs src)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT COMPONENTS RIO Tree TreePlayer MathCore Net Graf3d Graf Gpad EG REQUIRED)
find_package(LCIO)
find_package(podio)
find_package(plcio)
find_package(HepMC)
if(ROOT_FOUND)
    message("found ROOT: ${ROOT_INCLUDE_DIRS} ${ROOT_LIBRARIES}")
endif(ROOT_FOUND)
if(LCIO_FOUND)
    message("found LCIO: ${LCIO_INCLUDE_DIRS} ${LCIO_LIBRARIES}")
endif(LCIO_FOUND)
if(podio_FOUND)
    message("found podio: ${podio_INCLUDE_DIRS} ${podio_LIBRARIES}")
endif(podio_FOUND)
if(plcio_FOUND)
    message("found plcio: ${plcio_INCLUDE_DIRS} ${plcio_LIBRARY_DIR}")
endif(plcio_FOUND)
if(HepMC_FOUND)
    message("found HepMC: ${HepMC_INCLUDE_DIRS} ${HepMC_LIBRARY_DIR}")
endif(HepMC_FOUND)
############## for producing plcio library #############
INCLUDE_DIRECTORIES(${GenAlgo_incs})

gaudi_add_module(GenAlgo ${GenAlgo_srcs} 
  INCLUDE_DIRS 
    GaudiKernel
    FWCore
    Geant4
    ${LCIO_INCLUDE_DIRS}
    ${podio_INCLUDE_DIRS}
    ${plcio_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
    HepMC
  LINK_LIBRARIES 
    GaudiKernel 
    ${LCIO_LIBRARIES}
    ${podio_LIBRARIES}
    ROOT
    ${plcio_LIBRARY_DIR}/libplcio.so
    ${plcio_LIBRARY_DIR}/libplcioDict.so
    FWCore 
    HepMC
    Geant4
  )
#gaudi_add_test(Reader FRAMEWORK options/read.py)

###########################
