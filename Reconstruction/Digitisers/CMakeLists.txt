gaudi_subdir(Digitisers v0r0)

find_package(CLHEP REQUIRED)
find_package(GEAR REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(podio REQUIRED ) 
find_package(K4FWCore REQUIRED)

gaudi_depends_on_subdirs(
    Service/GearSvc
    Service/EventSeeder
)

# set(Digitisers_srcs
#     src/*.cpp
# )

# # Modules
# gaudi_add_module(Digitisers ${Digitisers_srcs}
#     INCLUDE_DIRS K4FWCore GaudiKernel GaudiAlgLib CLHEP gear ${plcio_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS}
#     LINK_LIBRARIES K4FWCore GaudiKernel GaudiAlgLib CLHEP $ENV{GEAR}/lib/libgearsurf.so ${GSL_LIBRARIES} $ENV{PLCIO}/lib/libplcio.so ${LCIO_LIBRARIES}
# )
