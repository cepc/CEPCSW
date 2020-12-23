
find_package(CLHEP REQUIRED;CONFIG)
find_package(EDM4HEP REQUIRED)
find_package(GSL REQUIRED )
message("GSL: ${GSL_LIBRARIES} ")
message("GSL INCLUDE_DIRS: ${GSL_INCLUDE_DIRS} ")

#gaudi_install_headers(DataHelper)

gaudi_add_library(DataHelperLib 
                  SOURCES src/CaloHitExtended.cc
                          src/Circle.cc
                          src/ClusterExtended.cc
                          src/ClusterShapes.cc
                          src/GroupTracks.cc
                          src/HelixClass.cc
                          src/LCCylinder.cc
                          src/LCLine3D.cc
                          src/LCPlane3D.cc
                          src/LineClass.cc
                          src/Navigation.cpp
                          src/SimpleHelix.cc
                          src/TrackerHitExtended.cc
                          src/TrackExtended.cc
                          src/TrackHitPair.cc
                  LINK EDM4HEP::edm4hep
                       EDM4HEP::edm4hepDict
                       ${GSL_LIBRARIES}
                       ${CLHEP_LIBRARIES}
)

install(TARGETS DataHelperLib
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
