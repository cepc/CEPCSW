gaudi_subdir(GearSvc v0r0)

find_package(GEAR REQUIRED)
find_package(DD4hep COMPONENTS DDCore DDRec REQUIRED)

gaudi_depends_on_subdirs(
    Detector/DetInterface
)

set(GearSvc_srcs
    src/*.cpp
)

gaudi_install_headers(GearSvc)

gaudi_add_module(GearSvc ${GearSvc_srcs}
    INCLUDE_DIRS GaudiKernel gear
    LINK_LIBRARIES GaudiKernel ${GEAR_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}
)
