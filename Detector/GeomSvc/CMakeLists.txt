
##############################################################################
# Package: GeomSvc
#    Desc: implement the IGeomSvc interface.
##############################################################################

gaudi_subdir(GeomSvc v0r0)

gaudi_depends_on_subdirs(
    FWCore
    Detector/DetInterface
)

find_package(ROOT COMPONENTS MathCore GenVector Geom Tree)
find_package(DD4hep COMPONENTS DDG4 DDRec REQUIRED)
find_package(GEAR REQUIRED)

message("GEAR_LIBRARIES: ${GEAR_LIBRARIES}")

gaudi_add_module(GeomSvc
                 src/GeomSvc.cpp
                 INCLUDE_DIRS
                   # DetInterface
                   # DD4hep
                   # GaudiKernel
                   # ROOT 
                 LINK_LIBRARIES
                   DD4hep 
                   ${DD4hep_COMPONENT_LIBRARIES} 
                   GaudiKernel
                   ${GEAR_LIBRARIES}
                   ${ROOT_LIBRARIES}
                   # ROOT
)
