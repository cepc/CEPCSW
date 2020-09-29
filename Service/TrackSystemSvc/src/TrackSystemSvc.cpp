#include "GearSvc/IGearSvc.h"
#include "gear/GearMgr.h"

#include "MarlinKalTest.h"

#include "TrackSystemSvc.h"

DECLARE_COMPONENT(TrackSystemSvc)

TrackSystemSvc::TrackSystemSvc(const std::string& name, ISvcLocator* svc)
  : base_class(name, svc),
    m_trackSystem(nullptr){
}

TrackSystemSvc::~TrackSystemSvc(){
}

MarlinTrk::IMarlinTrkSystem* TrackSystemSvc::getTrackSystem(){
  if(!m_trackSystem){
    gear::GearMgr* mgr=0; 
    auto _gear = service<IGearSvc>("GearSvc");
    if ( !_gear ) {
      info() << "Failed to find GearSvc ..." << endmsg;
    }
    else{
      mgr = _gear->getGearMgr();
    }

    auto _geoSvc = service<IGeoSvc>("GeoSvc");
    if ( !_geoSvc ) {
      info() << "Failed to find GeoSvc ..." << endmsg;
    }
    if(mgr==0&&_geoSvc==0){
      fatal() << "Both GearSvc and GeoSvc invalid!" << endmsg;
      return 0;
    }

    m_trackSystem = new MarlinTrk::MarlinKalTest( *mgr, _geoSvc ) ;
  }
  return m_trackSystem;
}

StatusCode TrackSystemSvc::initialize(){

  auto _gear = service<IGearSvc>("GearSvc");
  if ( !_gear ) {
    error() << "Failed to find GearSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  gear::GearMgr* mgr = _gear->getGearMgr();

  auto _geoSvc = service<IGeoSvc>("GeoSvc");
  if ( !_geoSvc ) {
    error() << "Failed to find GeoSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  m_trackSystem = new MarlinTrk::MarlinKalTest( *mgr, _geoSvc ) ;
  
  return StatusCode::SUCCESS;
}

void TrackSystemSvc::removeTrackSystem(){
  if ( m_trackSystem ) {
    delete m_trackSystem;
    m_trackSystem = nullptr;
  }
  return;
}

StatusCode TrackSystemSvc::finalize(){
  
  // if ( m_trackSystem ) {
  //   delete m_trackSystem;
  //   m_trackSystem = nullptr;
  // }
  
  return StatusCode::SUCCESS;
}
