/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "PlanarDigiAlg.h"
#include "DataHelper/TrackerHitHelper.h"
#include "GearSvc/IGearSvc.h"
#include "EventSeeder/IEventSeeder.h"
#include "TrackSystemSvc/ITrackSystemSvc.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/Vector3d.h"
/*
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include "UTIL/CellIDEncoder.h"
#include <UTIL/Operators.h>
*/
#include "Identifier/CEPCConf.h"
#include "UTIL/ILDConf.h"

// STUFF needed for GEAR
#include <gear/GEAR.h>
#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>

DECLARE_COMPONENT( PlanarDigiAlg )

PlanarDigiAlg::PlanarDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  //_description = "PlanarDigiAlg creates TrackerHits from SimTrackerHits, smearing them according to the input parameters." ;
  
  // Input collections
  declareProperty("HeaderCol", _headerCol);
  declareProperty("SimTrackHitCollection", _inColHdl, "Handle of the Input SimTrackerHit collection");
  
  // Output collections
  declareProperty("TrackerHitCollection", _outColHdl, "Handle of the TrackerHit output collection");
  declareProperty("TrackerHitAssociationCollection", _outRelColHdl, "Handle of TrackerHit SimTrackHit relation collection");
}

StatusCode PlanarDigiAlg::initialize()
{
  if (_parameterize) {
    if (_parU.size()!=10 || _parV.size()!=10) {
      fatal() << "parameters number must be 10! now " << _parU.size() << " for U and " << _parV.size() << " for V" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  else {
    if (_resU.size() !=  _resV.size()) {
      fatal() << "Inconsistent number of resolutions given for U and V coordinate: " 
              << "ResolutionU  :" <<   _resU.size() << " != ResolutionV : " <<  _resV.size()
              << endmsg;

      return StatusCode::FAILURE;
    }
  }

  // get the GEAR manager
  auto _gear = service<IGearSvc>("GearSvc");
  if ( !_gear ) {
    error() << "Failed to find GearSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  _GEAR = _gear->getGearMgr();

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  _SEEDER = service<IEventSeeder>("EventSeeder");
  _SEEDER->registerAlg(this);

  // zoujh: TODO - is this necessary? ...
  //FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
  //TODO:MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  //TODO:if( trksystem == 0 ) {
  //TODO:  throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
  //TODO:}
  //TODO:trksystem->init() ;
  //FIXME:SJA gear surface store has now been filled so we can dispose of the MarlinTrkSystem
  //TODO:delete trksystem;
  // fucd: fix TODO - surface store is needed to fill once always if does not handle tracking algorithm in job
  auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
  if ( !_trackSystemSvc ) {
    error() << "Failed to find TrackSystemSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  
  MarlinTrk::IMarlinTrkSystem* _trksystem =  _trackSystemSvc->getTrackSystem(this);
  _trksystem->init();
  
  _trackSystemSvc->removeTrackSystem(this);
  
  return GaudiAlgorithm::initialize();
}

StatusCode PlanarDigiAlg::execute()
{
  //auto header = _headerCol.get()->at(0);
  //int evtNo = header.getEventNumber();
  //int runNo = header.getRunNumber();
  //debug() << "Processing Run[" << runNo << "]::Event[" << evtNo << "]" << endmsg;

  //unsigned int thisSeed = _SEEDER->getSeed(this, evtNo, runNo);
  unsigned int thisSeed = _SEEDER->getSeed(this, _nEvt, 0);
  gsl_rng_set( _rng, thisSeed ) ;   
  debug() << "seed set to " << thisSeed << endmsg;

  auto trkhitVec = _outColHdl.createAndPut();
  auto relCol = _outRelColHdl.createAndPut();
  //auto STHcol = _inColHdl.get();
  //if ( STHcol == nullptr ) {
  //  debug() << "Collection " << _inColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  //  return StatusCode::SUCCESS;  // or FAILURE ?
  //}

  const edm4hep::SimTrackerHitCollection* STHcol = nullptr;
  try {
    STHcol = _inColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    return StatusCode::SUCCESS;
  }

  // **
  // ** STHcol != 0
  // **

  unsigned nCreatedHits=0;
  unsigned nDismissedHits=0;

  //CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::ILDCellID0::encoder_string , trkhitVec ) ;

  int nSimHits = STHcol->size()  ;

  int det_id = 0 ;
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;

  if ( nSimHits>0 ) {
    auto SimTHit = STHcol->at( 0 ) ;
    encoder.setValue(SimTHit.getCellID()) ;
    det_id  = encoder[lcio::ILDCellID0::subdet] ;
  
    if     ( det_id == lcio::ILDDetID::VXD ){}
    else if( det_id == lcio::ILDDetID::SIT ){}
    else if( det_id == lcio::ILDDetID::SET ){}
    else if( det_id == lcio::ILDDetID::FTD ){}
    else {
      fatal() << "unsupported detector ID = " << det_id << " CellID = " << SimTHit.getCellID()
              << ": file " << __FILE__ << " line " << __LINE__
              << endmsg;
      return StatusCode::FAILURE;
    }
  }

  //smearing
  debug() << " processing collection " << _inColHdl.fullKey()
          << " with " <<  nSimHits  << " hits ... " << endmsg ;

  int i = 0;
  for( auto SimTHit : *STHcol ) {
    if (SimTHit.getEDep()<=_eThreshold) continue;
    debug() << "MCParticle id " << SimTHit.getMCParticle().id() << endmsg;

    const int celId = SimTHit.getCellID() ;

    encoder.setValue(celId) ;
    int side   = encoder[lcio::ILDCellID0::side];
    int layer  = encoder[lcio::ILDCellID0::layer];
    int module = encoder[lcio::ILDCellID0::module];
    int sensor = encoder[lcio::ILDCellID0::sensor];

    debug() << "Hit = " << i << " has celId " << celId << endmsg;
    debug() << "side = " << side << endmsg;
    debug() << "layerNumber = " <<  layer << endmsg;
    debug() << "moduleNumber = " << module << endmsg;
    debug() << "sensorNumber = " << sensor << endmsg;

    //      //************************************************************
    //      // Quick check if the MCParticle is in the list of MCParticles
    //      //************************************************************
    //
    //      if ( ! SimTHit->getMCParticle().isAvailable() ) {
    //        error() << " SimHit " << SimTHit << " Created by zero MCParticle which is not in the list of MCParticles: "
    //                << endmsg;
    //        continue;
    //      }
    //
    //      //************************************************************
    //      // Quick check if the MCParticle is of zero charge
    //      //************************************************************
    //
    //      if( abs( SimTHit->getMCParticle().getCharge()) < 0.01  ){
    //        error() << " SimHit Created by zero charge particle: "
    //                << " Charge =  " << SimTHit->getMCParticle().getCharge()
    //                << " EDep =  " << SimTHit->getEDep()
    //                << " x =  " << SimTHit->getPosition()[0]
    //                << " PDG =  " << SimTHit->getMCParticle().getPDG()
    //                << " ID =  " << SimTHit->getMCParticle().id()
    //                << endmsg;
    //        continue;
    //      }
    gear::MeasurementSurface const* ms = _GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( encoder.lowWord() );
    gear::CartesianCoordinateSystem* cartesian = dynamic_cast< gear::CartesianCoordinateSystem* >( ms->getCoordinateSystem() );
    CLHEP::Hep3Vector uVec = cartesian->getLocalXAxis();
    CLHEP::Hep3Vector vVec = cartesian->getLocalYAxis();

    float u_direction[2] ;
    u_direction[0] = uVec.theta();
    u_direction[1] = uVec.phi();
    
    float v_direction[2] ;
    v_direction[0] = vVec.theta();
    v_direction[1] = vVec.phi();
    
    debug() << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1]
            << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
            << endmsg ;
    
    float resU(0), resV(0);

    if (!_parameterize) {
      if( (_resU.size() > 1 && layer > _resU.size()-1) || (_resV.size() > 1 && layer > _resV.size()-1) ) {
        fatal() << "layer exceeds resolution vector, please check input parameters ResolutionU and ResolutionV" << endmsg;
        return StatusCode::FAILURE;
      }

      resU = ( _resU.size() > 1 ?   _resU.value().at(  layer )     : _resU.value().at(0)   )  ;
      resV = ( _resV.size() > 1 ?   _resV.value().at(  layer )     : _resV.value().at(0)   )  ; 
    }
    else { // Riccardo's parameterized model
      auto mom = SimTHit.getMomentum();
      CLHEP::Hep3Vector momVec(mom[0], mom[1], mom[2]);
      const double alpha = uVec.azimAngle(momVec, vVec);
      const double cotanAlpha = 1./tan(alpha);
      // TODO: title angle (PI/2), magnetic field (3) 
      const double tanLorentzAngle = (side==0) ? 0. : 0.053 * 3 * cos(M_PI/2.);
      const double x = fabs(-cotanAlpha - tanLorentzAngle);
      resU = _parU[0] + _parU[1] * x + _parU[2] * exp(-_parU[9] * x) * cos(_parU[3] * x + _parU[4])
        + _parU[5] * exp(-0.5 * pow(((x - _parU[6]) / _parU[7]), 2)) + _parU[8] * pow(x, 0.5);
      
      const double beta = vVec.azimAngle(momVec, uVec);
      const double cotanBeta = 1./tan(beta);
      const double y = fabs(-cotanBeta);
      resV = _parV[0] + _parV[1] * y + _parV[2] * exp(-_parV[9] * y) * cos(_parV[3] * y + _parV[4])
        + _parV[5] * exp(-0.5 * pow(((y - _parV[6]) / _parV[7]), 2)) + _parV[8] * pow(y, 0.5);
    }

    debug() << " --- will smear hit with resU = " << resU << " and resV = " << resV << endmsg;

    auto& pos =  SimTHit.getPosition() ;  

    CLHEP::Hep3Vector globalPoint(pos[0],pos[1],pos[2]);
    CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
    CLHEP::Hep3Vector localPointSmeared = localPoint;

    debug() << std::setprecision(8) << "Position of hit before smearing global: ( "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " ) "
            << "local: ( " << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )" << endmsg;

    // A small check, if the hit is in the boundaries:
    if ( !ms->isLocalInBoundary( localPoint ) ){

      error() << "Hit local: ( " << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )"
              << " is not within boundaries. Hit is skipped." << endmsg;

      nDismissedHits++;
      continue;
    }

    unsigned  tries = 0;              

    bool accept_hit = false;

    // Now try to smear the hit and make sure it is still on the surface
    while( tries < 100 ) {

      if (tries > 0) {
        debug() << "retry smearing for side" << side << " layer"<< layer<< " module" << module
                << " sensor" << sensor << " : retries " << tries << endmsg;
      }

      double dx = gsl_ran_gaussian(_rng, resU);
      double dy = gsl_ran_gaussian(_rng, resV);
      localPointSmeared.setX( localPoint.x() + dx );
      localPointSmeared.setY( localPoint.y() + dy );

      //check if hit is in boundaries
      if ( ms->isLocalInBoundary( localPointSmeared ) && fabs(dx)<=_maxPull*resU && fabs(dy)<=_maxPull*resV ) {
        accept_hit = true;
        break;
      }

      tries++;
    }

    if( accept_hit == false ) {
      debug() << "hit could not be smeared within ladder after 100 tries: hit dropped"  << endmsg;
      continue; 
    }

    // for 1D strip measurements: set v to 0! Only the measurement in u counts!
    if( _isStrip || (resU!=0&&resV==0) ) localPointSmeared[1] = 0. ;

    // convert back to global position for TrackerHitPlaneImpl
    CLHEP::Hep3Vector globalPointSmeared = ms->getCoordinateSystem()->getGlobalPoint(localPointSmeared);

    debug() <<"Position of hit after smearing global: ( "  
            << globalPointSmeared.x() <<" "<< globalPointSmeared.y() <<" "<< globalPointSmeared.z() << " ) "
            << "local: ( "
            << localPointSmeared.x() << " " << localPointSmeared.y() << " " << localPointSmeared.z() << " )"
            << endmsg;
    
    //make the TrackerHitPlaneImpl
    auto trkHit = trkhitVec->create();

    trkHit.setCellID( encoder.lowWord() );

    edm4hep::Vector3d smearedPos(globalPointSmeared.x(), globalPointSmeared.y(), globalPointSmeared.z());
    trkHit.setPosition( smearedPos ) ;
    /*
    gear::CartesianCoordinateSystem* cartesian = dynamic_cast< gear::CartesianCoordinateSystem* >( ms->getCoordinateSystem() ); 
    CLHEP::Hep3Vector uVec = cartesian->getLocalXAxis();
    CLHEP::Hep3Vector vVec = cartesian->getLocalYAxis();

    float u_direction[2] ;
    u_direction[0] = uVec.theta();
    u_direction[1] = uVec.phi();

    float v_direction[2] ;
    v_direction[0] = vVec.theta();
    v_direction[1] = vVec.phi();

    debug() << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1]
            << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
            << endmsg ;
    */
    // fucd: next TODO: cov[0] = resU*reU, cov[2] = resV*resV, cov[5] = 0
    if(_usePlanarTag){
      std::array<float, 6> cov;
      cov[0] = u_direction[0];
      cov[1] = u_direction[1];
      cov[2] = resU;
      cov[3] = v_direction[0];
      cov[4] = v_direction[1];
      cov[5] = resV;
      trkHit.setCovMatrix(cov);
    /* zoujh: TODO - generate TrackerHitPlane with podio
    trkHit->setU( u_direction ) ;
    trkHit->setV( v_direction ) ;

    trkHit->setdU( resU ) ;

    if( _isStrip ) trkHit->setdV( 0 ); // no error in v direction for strip hits as there is no meesurement information in v direction
    else trkHit->setdV( resV ) ;
    */
      std::bitset<32> type;
      type.set(CEPCConf::TrkHitTypeBit::PLANAR);
      trkHit.setType((int)type.to_ulong());
    }
    else{
      trkHit.setCovMatrix(CEPC::ConvertToCovXYZ(resU, u_direction[0], u_direction[1], resV, v_direction[0], v_direction[1]));
    }

    if( _isStrip || (resU!=0&&resV==0) ){
        trkHit.setType( UTIL::set_bit( trkHit.getType() , UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) ) ;
    }
    trkHit.setEDep( SimTHit.getEDep() );
    trkHit.setTime( SimTHit.getTime() );
    // make the relation
    auto rel = relCol->create();

    float weight = 1.0;

    debug() <<" Set relation between "
            << " sim hit " << SimTHit.id() 
            << " to tracker hit " << trkHit.id()
            << " with a weight of " << weight 
            << endmsg;
    trkHit.addToRawHits(SimTHit.getObjectID());
    rel.setSim(SimTHit);
    rel.setRec(trkHit);
    rel.setWeight(weight);

    nCreatedHits++;

    debug() << "-------------------------------------------------------" << endmsg;

    ++i;
  }

  debug() << "Created " << nCreatedHits << " hits, "
          << nDismissedHits << " hits got dismissed for being out of boundary"
          << endmsg;
    
  _nEvt ++ ;

  return StatusCode::SUCCESS;
}

StatusCode PlanarDigiAlg::finalize()
{
  info() << "Processed " << _nEvt << " events " << endmsg;

  if(_rng) gsl_rng_free(_rng);

  return GaudiAlgorithm::finalize();
}
