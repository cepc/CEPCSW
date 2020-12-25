
gaudi_add_library(KiTrackLib 
                  SOURCES src/KiTrack/Automaton.cc 
                          src/KiTrack/HopfieldNeuralNet.cc 
                          src/KiTrack/IHit.cc 
                          src/KiTrack/SegmentBuilder.cc 
                          src/KiTrack/Segment.cc 

                          src/Criteria/Crit2_DeltaPhi.cc
                          src/Criteria/Crit2_DeltaPhi_MV.cc
                          src/Criteria/Crit2_DeltaRho.cc
                          src/Criteria/Crit2_DeltaTheta_MV.cc
                          src/Criteria/Crit2_Distance_MV.cc
                          src/Criteria/Crit2_HelixWithIP.cc
                          src/Criteria/Crit2_RZRatio.cc
                          src/Criteria/Crit2_StraightTrackRatio.cc
                          src/Criteria/Crit3_2DAngle.cc
                          src/Criteria/Crit3_2DAngleTimesR.cc
                          src/Criteria/Crit3_3DAngle.cc
                          src/Criteria/Crit3_3DAngleTimesR.cc
                          src/Criteria/Crit3_ChangeRZRatio.cc
                          src/Criteria/Crit3_IPCircleDist.cc
                          src/Criteria/Crit3_IPCircleDistTimesR.cc
                          src/Criteria/Crit3_NoZigZag_MV.cc
                          src/Criteria/Crit3_PT.cc
                          src/Criteria/Crit3_PT_MV.cc
                          src/Criteria/Crit4_2DAngleChange.cc
                          src/Criteria/Crit4_3DAngleChange.cc
                          src/Criteria/Crit4_3DAngleChangeNormed.cc
                          src/Criteria/Crit4_DistOfCircleCenters.cc
                          src/Criteria/Crit4_DistToExtrapolation.cc
                          src/Criteria/Crit4_NoZigZag.cc
                          src/Criteria/Crit4_PhiZRatioChange.cc
                          src/Criteria/Crit4_RChange.cc
                          src/Criteria/Criteria.cc
                          src/Criteria/SimpleCircle.cc

                          src/ILDImpl/FTDHit01.cc
                          src/ILDImpl/FTDHitSimple.cc
                          src/ILDImpl/FTDNeighborPetalSecCon.cc
                          src/ILDImpl/FTDSectorConnector.cc
                          src/ILDImpl/FTDTrack.cc
                          src/ILDImpl/MiniVector.cc
                          src/ILDImpl/MiniVectorHit01.cc
                          src/ILDImpl/SectorSystemFTD.cc
                          src/ILDImpl/SectorSystemVXD.cc
                          src/ILDImpl/VXDHit01.cc
                          src/ILDImpl/VXDHitSimple.cc
                          src/ILDImpl/VXDSectorConnector.cc
                          src/ILDImpl/VXDTrack.cc

                          src/Tools/Fitter.cc
                          src/Tools/FTDHelixFitter.cc
                          src/Tools/KiTrackMarlinTools.cc
                          src/Tools/Timer.cc

                  LINK DataHelperLib 
                       TrackSystemSvcLib 
                       ${ROOT_LIBRARIES}
                       ${CLHEP_LIBRARIES}
                       ${GSL_LIBRARIES}
                       EDM4HEP::edm4hep
                       ${LCIO_LIBRARIES}
)

target_include_directories(KiTrackLib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS KiTrackLib
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
