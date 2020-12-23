
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED )
find_package(EDM4HEP REQUIRED ) 
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)
find_package(k4FWCore)

# Modules
gaudi_add_module(SiliconTracking
                 SOURCES src/ForwardTrackingAlg.cpp
                         src/SiliconTrackingAlg.cpp
                         src/SpacePointBuilderAlg.cpp
                         src/TrackSubsetAlg.cpp
                 LINK GearSvc
                      EventSeeder
                      TrackSystemSvcLib
                      DataHelperLib
                      KiTrackLib
                      Gaudi::GaudiKernel
                      k4FWCore::k4FWCore 
                      ${GEAR_LIBRARIES} 
                      ${GSL_LIBRARIES} 
                      ${LCIO_LIBRARIES} 
)
install(TARGETS SiliconTracking
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
