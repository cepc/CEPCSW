
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
			  src/TrackerHitHelper.cpp
                  LINK EDM4HEP::edm4hep
                       EDM4HEP::edm4hepDict
                       ${GSL_LIBRARIES}
                       ${CLHEP_LIBRARIES}
		       Identifier
)

install(TARGETS DataHelperLib
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
