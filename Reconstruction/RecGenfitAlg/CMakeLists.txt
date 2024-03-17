# RecGenfitAlg

if (GenFit_FOUND)
gaudi_add_module(RecGenfitAlg
       SOURCES src/RecGenfitAlgSDT.cpp
                src/GenfitTrack.cpp
                src/GenfitHit.cpp
                src/GenfitField.cpp
                src/GenfitFitter.cpp
                src/GenfitMaterialInterface.cpp
                src/LSFitting.cpp
                src/WireMeasurementDC.cpp
                src/PlanarMeasurementSDT.cpp
        LINK GearSvc
             Gaudi::GaudiAlgLib
             Gaudi::GaudiKernel
             ${GEAR_LIBRARIES}
             ${GSL_LIBRARIES}
             ${LCIO_LIBRARIES}
             DetSegmentation
             DetInterface
             DataHelperLib
             EDM4HEP::edm4hep
             EDM4HEP::edm4hepDict
             GenFit::genfit2
             k4FWCore::k4FWCore
)

target_include_directories(RecGenfitAlg PUBLIC
        ${LCIO_INCLUDE_DIRS}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS RecGenfitAlg
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
endif()
