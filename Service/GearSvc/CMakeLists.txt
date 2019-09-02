gaudi_subdir(GearSvc v0r0)

find_package(GEAR REQUIRED)

set(GearSvc_srcs
    src/*.cpp
)

gaudi_install_headers(GearSvc)

gaudi_add_module(GearSvc ${GearSvc_srcs}
    INCLUDE_DIRS GaudiKernel gear
    LINK_LIBRARIES GaudiKernel $ENV{GEAR}/lib/libgear.so $ENV{GEAR}/lib/libgearxml.so
)
