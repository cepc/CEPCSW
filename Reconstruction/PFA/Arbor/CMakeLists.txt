
# Modules
gaudi_add_module(Arbor
                 SOURCES src/Arbor.cc
		 	src/ArborHit.cc
		 	src/ArborTool.cc
		 	src/ArborToolLCIO.cc 
		 	src/ClusterAna.cc 
		 	src/DetectorPos.cc 
		 	src/HelixClassD.cc 
			 src/MarlinArbor.cc
		 	src/BushConnect.cc
                 LINK EventSeeder
                      GearSvc
                      DataHelperLib
                      DetInterface
                      Gaudi::GaudiKernel
                      k4FWCore::k4FWCore
                      ${LCContent_LIBRARIES}
                      ${CLHEP_LIBRARIES} 
                      ${ROOT_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      ${GEAR_LIBRARIES}
                      ${DD4hep_COMPONENT_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)


install(TARGETS Arbor
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)

