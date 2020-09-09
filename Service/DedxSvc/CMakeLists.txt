gaudi_subdir(DedxSvc v0r0)

set(DedxSvc_srcs
    src/*.cpp
)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
gaudi_install_headers(DedxSvc)

gaudi_add_module(DedxSvc ${DedxSvc_srcs}
    INCLUDE_DIRS GaudiKernel
    LINK_LIBRARIES GaudiKernel
)
