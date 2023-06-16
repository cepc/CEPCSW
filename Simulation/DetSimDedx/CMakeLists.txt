
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
find_package(Garfield REQUIRED)
message("libonnxruntime ${OnnxRuntime_LIBRARY}")
find_package(OnnxRuntime REQUIRED)

message("libonnxruntime ${OnnxRuntime_LIBRARY}")
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
                      k4FWCore::k4FWCore
                      Garfield::Garfield
                      ${OnnxRuntime_LIBRARY}
                      #/cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/lib/libonnxruntime.so
                      ${CLHEP_LIBRARIES}
 
)

