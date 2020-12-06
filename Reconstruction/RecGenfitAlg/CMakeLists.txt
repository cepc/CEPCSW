gaudi_subdir(RecGenfitAlg v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(GSL REQUIRED )
find_package(LCIO REQUIRED )
find_package(podio REQUIRED )
find_package(EDM4HEP REQUIRED)
find_package(ROOT REQUIRED Geom)
find_package(DD4hep COMPONENTS DDCore DDRec DDParsers REQUIRED)

gaudi_depends_on_subdirs(
    Detector/DetInterface
    Detector/DetSegmentation
    Utilities/DataHelper
)

set (GenFit_INCLUDE_DIRS $ENV{GENFIT_ROOT}/include)
set (GenFit_LIB_DIRS $ENV{GENFIT_ROOT}/lib64)
set (Eigen_INCLUDE_DIRS $ENV{Eigen_ROOT}/include/eigen3)

set(RecGenfitAlg_srcs
  src/RecGenfitAlgDC.cpp
  src/GenfitTrack.cpp
  src/GenfitField.cpp
  src/GenfitFitter.cpp
  src/GenfitMaterialInterface.cpp
  src/GenfitMsg.cpp
  )

# Modules
gaudi_add_module(RecGenfitAlg ${RecGenfitAlg_srcs}
  INCLUDE_DIRS k4FWCore GaudiKernel GaudiAlgLib CLHEP ROOT gear
  ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS} ${GenFit_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  LINK_LIBRARIES k4FWCore GaudiKernel GaudiAlgLib CLHEP ROOT DataHelperLib
  DetSegmentation $ENV{GEAR}/lib/libgearsurf.so ${GSL_LIBRARIES} ${LCIO_LIBRARIES}
  ${GenFit_LIB_DIRS}/libgenfit2.so DD4hep ${DD4hep_COMPONENT_LIBRARIES}
  -Wl,--no-as-needed EDM4HEP::edm4hep EDM4HEP::edm4hepDict -Wl,--as-needed
  )
