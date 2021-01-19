# RecGenfitAlg

if (GenFit_FOUND)
gaudi_add_module(RecGenfitAlg
        SOURCES src/RecGenfitAlgDC.cpp
                src/GenfitTrack.cpp
                src/GenfitField.cpp
                src/GenfitFitter.cpp
                src/GenfitMaterialInterface.cpp
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
)

target_include_directories(RecGenfitAlg PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS RecGenfitAlg
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
endif()
