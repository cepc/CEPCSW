gaudi_subdir(Event v0r0)

# Declare PODIO dependency
find_package(podio REQUIRED HINTS $ENV{PODIO})
## Make sure the library is found.
## Not the case if LD_LIBRARY_PATH is wrong
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

# Declare ROOT dependency
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(ROOT REQUIRED COMPONENTS RIO Tree)
include_directories(${ROOT_INCLUDE_DIR})
include(${ROOT_USE_FILE})

# Settings
file(GLOB sources src/*.cc)
file(GLOB headers Event/*.h)
include_directories(Event)

# Dictionary and Modules
REFLEX_GENERATE_DICTIONARY(Event ${headers} SELECTION src/selection.xml)
add_library(EventDict SHARED Event.cxx)
add_dependencies(EventDict Event-dictgen)
target_link_libraries(EventDict Event ROOT::Core)

add_library(Event SHARED ${sources} ${headers})
target_link_libraries(Event podio ROOT::RIO ROOT::Tree)

# Headers and Libraries
gaudi_install_headers(Event)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/EventDict.rootmap DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Event_rdict.pcm DESTINATION lib)

# Unit tests
#gaudi_add_test(HelloAlg
#               FRAMEWORK options/helloalg.py)
