# RecGenfitAlg
gaudi_add_module(RecGenfitAlg
        SOURCES src/RecGenfitAlgSDT.cpp
                src/RecGenfitAlgDC.cpp
                src/GenfitTrack.cpp
                src/GenfitField.cpp
                src/GenfitFitter.cpp
                src/GenfitMaterialInterface.cpp
                src/GenfitMsg.cpp
        LINK GearSvc
#             EventSeeder
#             TrackSystemSvcLib
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
             $ENV{GENFIT_ROOT}/lib64/libgenfit2.so

)

target_include_directories(RecGenfitAlg PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:$ENV{GENFIT_ROOT}>/include
)

install(TARGETS RecGenfitAlg
            EXPORT CEPCSWTargets
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
            COMPONENT dev)
