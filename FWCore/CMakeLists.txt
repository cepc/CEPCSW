################################################################################
# Package: FWCore
################################################################################
gaudi_subdir(FWCore v0r1)

find_package(podio REQUIRED)
find_package(plcio REQUIRED)
find_package(LCIO REQUIRED)

# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg GaudiKernel)
# gaudi_install_scripts()
gaudi_install_python_modules()

gaudi_add_library(FWCore
		  src/*.cpp
                  INCLUDE_DIRS podio ${LCIO_INCLUDE_DIRS} ${plcio_INCLUDE_DIRS} ROOT
                  LINK_LIBRARIES GaudiAlgLib GaudiKernel podio::podioRootIO ${LCIO_LIBRARIES} ${plcio_LIBRARIES} ROOT
                  PUBLIC_HEADERS FWCore)

gaudi_add_module(FWCorePlugins
                 src/components/*.cpp
                 LINK_LIBRARIES GaudiAlgLib GaudiKernel FWCore ROOT)

