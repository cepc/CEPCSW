
find_package(CLHEP REQUIRED;CONFIG)
find_package(EDM4HEP REQUIRED)
find_package(GSL REQUIRED ) 
find_package(LCIO REQUIRED ) 
find_package(k4FWCore REQUIRED)

set(TotalInvMass_srcs src/TotalInvMass.cc)

gaudi_add_module(TotalInvMass
                 SOURCES ${TotalInvMass_srcs}
                 LINK k4FWCore::k4FWCore 
                      Gaudi::GaudiKernel 
                      Gaudi::GaudiAlgLib 
                      ${CLHEP_LIBRARIES} 
                      ${GSL_LIBRARIES}
                      ${LCIO_LIBRARIES}
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict
                      ${ROOT_LIBRARIES}
)
