gaudi_subdir(Event v0r0)

find_package(podio REQUIRED HINTS $ENV{PODIO})
# Make sure the library is found.
# Not the case if LD_LIBRARY_PATH is wrong
find_library(PODIOLIB podio PATHS ${podio_DIR}/../lib)
if (NOT PODIOLIB)
    message(FATAL_ERROR "libpodio.so(dylib) cannot be found dynamically. Make sure you have set up your environment to use PODIO")
endif()
link_directories(${podio_LIBRARY_DIR})

execute_process(COMMAND python ${podio_CMAKE_DIR}/../python/podio_class_generator.py Event/plcio.yaml Event Event
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

include_directories(
    ${podio_INCLUDE_DIRS}
    )

set(Event_srcs
    src/*.cc
    )

# Headers and Libraries
gaudi_install_headers(Event)

# Modules
gaudi_add_module(Event ${Event_srcs}
    INCLUDE_DIRS Event
    #LINK_LIBRARIES GaudiKernel
)

# Unit tests
#gaudi_add_test(HelloAlg
#               FRAMEWORK options/helloalg.py)
