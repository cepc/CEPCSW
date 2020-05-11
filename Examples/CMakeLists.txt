
gaudi_subdir(Examples v0r0)

find_package(podio REQUIRED)
find_package(plcio REQUIRED)
find_package(LCIO REQUIRED)
find_package(EDM4HEP REQUIRED)

set(Examples_srcs
    src/HelloWorld/*.cpp
    src/FirstSvc/*.cpp
    src/SecondAlg/*.cpp
    src/PlcioTest/*.cpp
    src/Edm4hepTest/*.cpp
)

# Headers and Libraries

gaudi_install_headers(Examples)

# Modules
gaudi_add_module(Examples ${Examples_srcs}
    INCLUDE_DIRS GaudiKernel FWCore ${plcio_INCLUDE_DIRS} ${podio_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
    LINK_LIBRARIES GaudiKernel FWCore ${podio_LIBRARIES} ${LCIO_LIBRARIES} $ENV{PLCIO}/lib/libplcio.so 
      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)

# Unit tests
gaudi_add_test(HelloAlg
               FRAMEWORK options/helloalg.py)

gaudi_add_test(SecondAlg
               FRAMEWORK options/secondalg.py)

gaudi_add_test(PlcioWriteAlg
               FRAMEWORK options/plcio_write.py)

gaudi_add_test(PlcioReadAlg
               FRAMEWORK options/plcio_read.py)

gaudi_add_test(Edm4hepWriteAlg
               FRAMEWORK options/edm4hep_write.py)

gaudi_add_test(Edm4hepReadAlg
               FRAMEWORK options/edm4hep_read.py)

gaudi_add_test(LCIOReadAlg
               FRAMEWORK options/LCIO_read.py)

