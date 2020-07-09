gaudi_subdir(GaudiPandora v0r0)

find_package(CLHEP REQUIRED)
#find_package(GSL REQUIRED )##don't use this,  use CEPC LCIO version one , due to the ClusterShape.cc file which is from LCIO
#message("GSL: ${GSL_LIBRARIES} ")
set (gsl_include "/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/GSL/1.14/install/include")
set (gsl_lib1 "/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/GSL/1.14/install/lib/libgsl.so")
set (gsl_lib2 "/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/GSL/1.14/install/lib/libgslcblas.so")
find_package(LCIO REQUIRED ) 
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


list(APPEND CMAKE_MODULE_PATH "$ENV{ROOTSYS}/etc/cmake/")
find_package(ROOT 5.26.00 REQUIRED COMPONENTS Eve Geom RGL EG)

#In order to use ClusterShape.cc
include_directories("../CED/CED/")
include_directories("../MarlinUtil/01-08/source/")

gaudi_depends_on_subdirs(
    Service/EventSeeder
    Service/GearSvc
)

set(dir_srcs
    src/PandoraPFAlg.cpp
    src/MCParticleCreator.cpp
    src/GeometryCreator.cpp
    src/CaloHitCreator.cpp
    src/TrackCreator.cpp
    src/PfoCreator.cpp
    ../CED/CED/*.cc
    ../MarlinUtil/01-08/source/*.cc
)
set(dir_include include)
# Modules
gaudi_add_module(GaudiPandora ${dir_srcs}
    INCLUDE_DIRS ${gsl_include} ${dir_include} GaudiKernel FWCore CLHEP  ${LCIO_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} gear  
    LINK_LIBRARIES ${gsl_lib1} ${gsl_lib2} GaudiKernel FWCore CLHEP ROOT ${LCIO_LIBRARIES} $ENV{GEAR}/lib/libgear.so $ENV{GEAR}/lib/libgearxml.so 
      -Wl,--no-as-needed 
      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
      -Wl,--as-needed 

)
