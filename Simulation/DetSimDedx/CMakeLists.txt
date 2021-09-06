
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
#include_directories("/cvmfs/cepcsw.ihep.ac.cn/prototype/releases/externals/98.0.0/Garfield++/include/")
find_package(Garfield REQUIRED)
gaudi_add_module(DetSimDedx
                 SOURCES src/DummyDedxSimTool.cpp
                         src/BetheBlochEquationDedxSimTool.cpp
                         src/GFDndxSimTool.cpp
                         src/TrackHeedSimTool.cpp

                 LINK DetSimInterface
                      DetInterface
                      DetSegmentation
                      ${DD4hep_COMPONENT_LIBRARIES}
                      Gaudi::GaudiKernel
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      Garfield::Garfield
                      #/cvmfs/cepcsw.ihep.ac.cn/prototype/releases/externals/98.0.0/Garfield++/lib64/libGarfield.so
)
