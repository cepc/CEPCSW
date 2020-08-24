gaudi_subdir(KiTrack v0r0)

find_package(ROOT REQUIRED)
find_package(CLHEP REQUIRED)
#find_package(DD4hep REQUIRED)
find_package(GSL REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(LCIO REQUIRED)

gaudi_depends_on_subdirs(Service/TrackSystemSvc Utilities/DataHelper)

set(KiTrackLib_srcs src/KiTrack/*.cc src/Criteria/*.cc src/ILDImpl/*.cc src/Tools/*.cc)

#gaudi_install_headers(src)
include_directories(src)

gaudi_add_library(KiTrackLib ${KiTrackLib_srcs}
    PUBLIC_HEADERS KiTrack 
    LINK_LIBRARIES DataHelperLib TrackSystemSvcLib ROOT CLHEP GSL EDM4HEP::edm4hep LCIO
    #              DD4hep
)
