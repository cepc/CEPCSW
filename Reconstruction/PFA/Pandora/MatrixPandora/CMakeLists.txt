gaudi_subdir(MatrixPandora v0r0)

find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(CLHEP REQUIRED)
find_package(LCIO REQUIRED ) 
find_package(GEAR REQUIRED)
find_package(EDM4HEP REQUIRED ) 
include_directories(${EDM4HEP_INCLUDE_DIR})

find_package(PandoraSDK REQUIRED ) 
find_package(LCContent REQUIRED ) 
include_directories(${PandoraSDK_INCLUDE_DIRS})
link_libraries(${PandoraSDK_LIBRARIES})
include_directories(${LCContent_INCLUDE_DIRS})
link_libraries(${LCContent_LIBRARIES})


list(APPEND CMAKE_MODULE_PATH "$ENV{ROOTSYS}/etc/cmake/")
find_package(ROOT 5.26.00 REQUIRED COMPONENTS Eve Geom RGL EG)


gaudi_depends_on_subdirs(
    Service/EventSeeder
    Service/GearSvc
    Detector/DetInterface
    Utilities/DataHelper
)

set(dir_srcs
    src/PandoraMatrixAlg.cpp
    src/MCParticleCreator.cpp
    src/GeometryCreator.cpp
    src/CaloHitCreator.cpp
    src/TrackCreator.cpp
    src/PfoCreator.cpp
)
set(dir_include include)
# Modules
gaudi_add_module(MatrixPandora ${dir_srcs}
    INCLUDE_DIRS ${dir_include} GaudiKernel FWCore CLHEP  ${LCIO_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} gear DD4hep  
    LINK_LIBRARIES GaudiKernel FWCore CLHEP ROOT ${LCIO_LIBRARIES} $ENV{GEAR}/lib/libgear.so $ENV{GEAR}/lib/libgearxml.so DD4hep ${DD4hep_COMPONENT_LIBRARIES} DDRec DataHelperLib
      -Wl,--no-as-needed 
      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
      -Wl,--as-needed 
)
