gaudi_subdir(DataHelper v0r0)

find_package(CLHEP REQUIRED;CONFIG)
find_package(EDM4HEP REQUIRED)
find_package(GSL REQUIRED )
message("GSL: ${GSL_LIBRARIES} ")
message("GSL INCLUDE_DIRS: ${GSL_INCLUDE_DIRS} ")

# gaudi_depends_on_subdirs()

set(DataHelperLib_srcs src/*.cc src/*.cpp)

#gaudi_install_headers(DataHelper)

gaudi_add_library(DataHelperLib ${DataHelperLib_srcs}
    PUBLIC_HEADERS DataHelper
    INCLUDE_DIRS ${CLHEP_INCLUDE_DIR} ${GSL_INCLUDE_DIRS}
    LINK_LIBRARIES EDM4HEP::edm4hep EDM4HEP::edm4hepDict ${GSL_LIBRARIES} ${CLHEP_LIBRARIES}
)
