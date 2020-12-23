gaudi_subdir(TrackInspect v0r0)

find_package(EDM4HEP REQUIRED) 
find_package(AIDA REQUIRED)
find_package(ROOT REQUIRED)
find_package(CLHEP REQUIRED;CONFIG)

gaudi_depends_on_subdirs(
    Utilities/DataHelper
)

set(TrackInspect_srcs src/*.cpp)

# Modules
gaudi_add_module(TrackInspect ${TrackInspect_srcs}
    INCLUDE_DIRS AIDA ROOT
    LINK_LIBRARIES DataHelperLib GaudiKernel EDM4HEP ROOT ${CLHEP_LIBRARIES}
)
