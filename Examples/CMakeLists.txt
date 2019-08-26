
gaudi_subdir(Examples v0r0)

find_package(podio REQUIRED)
find_package(plcio REQUIRED)

set(Examples_srcs
    src/HelloWorld/*.cpp
    src/FirstSvc/*.cpp
    src/SecondAlg/*.cpp
    src/PlcioTest/*.cpp
)

# Headers and Libraries

gaudi_install_headers(Examples)

# Modules
gaudi_add_module(Examples ${Examples_srcs}
    INCLUDE_DIRS GaudiKernel FWCore ${plcio_INCLUDE_DIRS} ${podio_INCLUDE_DIRS}
    LINK_LIBRARIES GaudiKernel FWCore ${podio_LIBRARIES} $ENV{PLCIO}/lib/libplcio.so
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
