gaudi_add_module(DumpEvent
                 SOURCES src/DumpMCParticleAlg.cpp
                         src/DumpSimHitAlg.cpp
		         #src/DumpHitAlg.cpp
			 src/DumpTrackAlg.cpp
			 #src/DumpCalorimeterAlg.cpp
		 LINK DataHelperLib
                      Gaudi::GaudiKernel
                      EDM4HEP::edm4hep
                      ${ROOT_LIBRARIES}
                      ${CLHEP_LIBRARIES}
		      ${DD4hep_COMPONENT_LIBRARIES}
                      DetInterface
                      k4FWCore::k4FWCore
)

install(TARGETS DumpEvent
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
