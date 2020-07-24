gaudi_subdir(DataHelper v0r0)

find_package(EDM4HEP REQUIRED)

gaudi_depends_on_subdirs()

set(DataHelperLib_srcs src/*.cc src/*.cpp)

#gaudi_install_headers(DataHelper)

gaudi_add_library(DataHelperLib ${DataHelperLib_srcs}
    PUBLIC_HEADERS DataHelper
    LINK_LIBRARIES EDM4HEP::edm4hep EDM4HEP::edm4hepDict
)
