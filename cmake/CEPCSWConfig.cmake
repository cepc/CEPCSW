include(CMakeFindDependencyMacro)
find_dependency(podio REQUIRED)
find_dependency(Gaudi REQUIRED)
find_dependency(k4FWCore REQUIRED)
find_dependency(k4LCIOReader REQUIRED)
find_dependency(EDM4HEP REQUIRED)
find_dependency(ROOT REQUIRED)

# - Include the targets file to create the imported targets that a client can
# link to (libraries) or execute (programs)
include("${CMAKE_CURRENT_LIST_DIR}/CEPCSWTargets.cmake")

get_property(TEST_CEPCSW_LIBRARY TARGET CEPCSW::GeomSvc PROPERTY LOCATION)
find_package_handle_standard_args(CEPCSW DEFAULT_MSG CMAKE_CURRENT_LIST_FILE TEST_CEPCSW_LIBRARY)
