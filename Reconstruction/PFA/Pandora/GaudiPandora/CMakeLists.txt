gaudi_subdir(GaudiPandora v0r0)

find_package(LCIO REQUIRED ) 
find_package(DD4hep COMPONENTS DDG4 REQUIRED)
find_package(GEAR REQUIRED)
message("ENV GEAR: $ENV{GEAR}")

find_package(EDM4HEP REQUIRED )
include_directories(${EDM4HEP_INCLUDE_DIR})

find_package(PandoraSDK REQUIRED ) 
include_directories(${PandoraSDK_INCLUDE_DIRS})
link_libraries(${PandoraSDK_LIBRARIES})
find_package(LCContent REQUIRED ) 
include_directories(${LCContent_INCLUDE_DIRS})
link_libraries(${LCContent_LIBRARIES})


find_package(ROOT COMPONENTS MathCore Physics GenVector Geom REQUIRED)


gaudi_depends_on_subdirs(
    Service/EventSeeder
    Service/GearSvc
    Utilities/DataHelper
    Detector/DetInterface
)

set(dir_srcs
    src/PandoraPFAlg.cpp
    src/MCParticleCreator.cpp
    src/GeometryCreator.cpp
    src/CaloHitCreator.cpp
    src/TrackCreator.cpp
    src/PfoCreator.cpp
    src/Utility.cpp
)
set(dir_include include)
# Modules
gaudi_add_module(GaudiPandora ${dir_srcs}
    INCLUDE_DIRS ${dir_include} GaudiKernel FWCore CLHEP  ${LCIO_INCLUDE_DIRS} ROOT gear DD4hep 
    LINK_LIBRARIES GaudiKernel FWCore CLHEP ROOT ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} DataHelperLib DD4hep ${DD4hep_COMPONENT_LIBRARIES} DDRec 
      -Wl,--no-as-needed 
      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
      -Wl,--as-needed 

)
