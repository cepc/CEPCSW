#include "SpacePointBuilderAlg.h"
#include "GearSvc/IGearSvc.h"
#include "TrackSystemSvc/ITrackSystemSvc.h"

#include "UTIL/ILDConf.h"
//#include "UTIL/LCRelationNavigator.h"

#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"
#include "gear/ZPlanarParameters.h"
#include "gear/ZPlanarLayerLayout.h"

#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"

#include <cmath>
#include <sstream>

DECLARE_COMPONENT(SpacePointBuilderAlg)

SpacePointBuilderAlg::SpacePointBuilderAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc) {
  //declareProperty("HeaderCol", _headerColHdl);
  declareProperty("MCParticleCollection", _inMCColHdl, "Handle of the Input MCParticle collection");
  declareProperty("TrackerHitCollection", _inHitColHdl, "Handle of the Input TrackerHits collection");
  //declareProperty("FTDPixelHitCollection", _inFTDPixelColHdl, "Handle of the Input FTD TrackerHits collection");
  declareProperty("TrackerHitAssociationCollection", _inHitAssColHdl, "Handle of the Input MCRecoTrackerAssociation collection");
  //declareProperty("SITHitCollection", _inSITColHdl, "Handle of the Input SIT TrackerHits collection");
  declareProperty("SpacePointCollection", _outSPColHdl, "Handle of the SpacePoint output collection");
  declareProperty("SpacePointAssociationCollection", _outSPAssColHdl, "Handle of the SpacePoints association output collection");
}

StatusCode SpacePointBuilderAlg::initialize() { 
  debug() << "   init called  " << endmsg;
  // usually a good idea to
  _nRun = 0 ;
  _nEvt = 0 ;
  
  _nominal_vertex.set(_nominal_vertex_x, _nominal_vertex_y, _nominal_vertex_z);

  auto gearSvc = service<IGearSvc>("GearSvc");
  if ( !gearSvc ) {
    error() << "Failed to find GearSvc ..." << endmsg;
    return StatusCode::FAILURE;
  }
  _GEAR = gearSvc->getGearMgr();
  
  //FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
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

StatusCode SpacePointBuilderAlg::execute(){
  StatusCode sc;
  auto spCol = _outSPColHdl.createAndPut();
  auto relCol = _outSPAssColHdl.createAndPut();

  const edm4hep::TrackerHitCollection* hitCol = nullptr;
  try {
    hitCol = _inHitColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inHitColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    return sc;
  }

  const edm4hep::MCRecoTrackerAssociationCollection* hitAssCol = nullptr;
  try {
    hitAssCol = _inHitAssColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inHitAssColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    return sc;
  }
    
  if(hitCol!=NULL && hitAssCol!=NULL){
    unsigned createdSpacePoints = 0;
    unsigned rawStripHits = 0;
    unsigned possibleSpacePoints = 0;
    _nOutOfBoundary = 0;
    _nStripsTooParallel = 0;
    _nPlanesNotParallel = 0;
    
    //edm4hep::TrackerHitCollection* spCol = new edm4hep::TrackerHitCollection();    // output spacepoint collection
    //edm4hep::MCRecoTrackerAssociationCollection* relCol = new edm4hep::MCRecoTrackerAssociationCollection();    // output relation collection
    
    // to store the weights
    //LCFlagImpl lcFlag(0) ;
    //lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    //relCol->setFlag( lcFlag.getFlag()  ) ;
    
    unsigned nHits = hitCol->size();
    debug() << "Number of hits: " << nHits << endmsg;
    
    //store hits in map according to their CellID0
    std::map<long long, std::vector<edm4hep::TrackerHit> > map_cellID0_hits;
    std::map<long long, std::vector<edm4hep::TrackerHit> >::iterator it;
    for(auto trkHit : *hitCol){
      debug() << "Add hit with CellID0 = " << trkHit.getCellID() << " " << getCellID0Info( trkHit.getCellID() ) << endmsg;
      map_cellID0_hits[ trkHit.getCellID() ].push_back( trkHit );
    }
    
    UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
    // now loop over all CellID0s
    for( it= map_cellID0_hits.begin(); it!= map_cellID0_hits.end(); it++ ){
      rawStripHits += it->second.size();
      std::vector<edm4hep::TrackerHit>& hitsFront = it->second;
      unsigned long long cellID0 = it->first;
      //get the CellID0s at the back of this sensor
      std::vector<int> cellID0sBack = getCellID0sAtBack( cellID0 );
      
      for( unsigned i=0; i< cellID0sBack.size(); i++ ){ 
        int cellID0Back = cellID0sBack[i];
        std::vector<edm4hep::TrackerHit>& hitsBack = map_cellID0_hits[ cellID0Back ];
        debug() << "strips: CellID0 " << cellID0  << " " << getCellID0Info( cellID0 )  << "(" << hitsFront.size()
		<< " hits) <---> CellID0 " << cellID0Back << getCellID0Info( cellID0Back )
		<< "(" << hitsBack.size() << " hits)" << endmsg
		<< "--> " << hitsFront.size() * hitsBack.size() << " possible combinations" << endmsg;
        
        possibleSpacePoints += hitsFront.size() * hitsBack.size();
	// Now iterate over all combinations and store those that make sense
        for( unsigned ifront=0; ifront<hitsFront.size(); ifront++ ){
	  edm4hep::TrackerHit& hitFront = hitsFront[ifront];
	  for( unsigned j=0; j<hitsBack.size(); j++ ){
	    edm4hep::TrackerHit& hitBack = hitsBack[j];

	    std::vector<edm4hep::ConstSimTrackerHit> simHitsFront;
	    std::vector<edm4hep::ConstSimTrackerHit> simHitsBack;
	    for(auto hitAss : *hitAssCol){
	      if(hitAss.getRec().id()==hitFront.id()) simHitsFront.push_back(hitAss.getSim());
	      if(hitAss.getRec().id()==hitBack.id()) simHitsBack.push_back(hitAss.getSim());
	    }
	    debug() << "attempt to create space point from:" << endmsg;
            debug() << "   front hit: " << hitFront.id() << " no. of simhit = " << simHitsFront.size() ;
            if( simHitsFront.size()!=0 ) { 
	      edm4hep::ConstSimTrackerHit& simhit = simHitsFront[0];
              debug() << "   first simhit = " << simhit.id() << " mcp = " << simhit.getMCParticle().id() << " (" << simhit.getPosition() << ") " ; 
            }
            debug() << endmsg;
            debug() << "  rear hit: " << hitBack.id() << " no. of simhit = " << simHitsBack.size() ;
            if( simHitsBack.size()!=0 ) { 
	      edm4hep::ConstSimTrackerHit& simhit = simHitsBack[0];
              debug() << "   first simhit = " << simhit.id() << " mcp = "<< simhit.getMCParticle().id() << " (" << simhit.getPosition() << ") " ; 
            }
	    debug() << endmsg;
            
            bool ghost_hit = true;
	    if (simHitsFront.size()==1 && simHitsBack.size() == 1) {
              debug() << "SpacePoint creation from two good hits:" << endmsg;
	      ghost_hit = simHitsFront[0].getMCParticle().id() != simHitsBack[0].getMCParticle().id();
	    }
	    if ( ghost_hit == true ) {
              debug() << "SpacePoint Ghosthit!" << endmsg;
            }
            
            cellID.setValue( cellID0 );
	    int subdet = cellID[ UTIL::ILDCellID0::subdet ] ;
            double strip_length_mm = 0;
            if (subdet == UTIL::ILDDetID::SIT) {
	      strip_length_mm = _GEAR->getSITParameters().getDoubleVal("strip_length_mm");
            }
	    else if (subdet == UTIL::ILDDetID::SET) {
              strip_length_mm = _GEAR->getSETParameters().getDoubleVal("strip_length_mm");
            }
	    else if (subdet == UTIL::ILDDetID::FTD) {
              strip_length_mm = _GEAR->getFTDParameters().getDoubleVal("strip_length_mm");
            }
	    else {
	      std::stringstream errorMsg;
              errorMsg << "SpacePointBuilderAlg::processEvent: unsupported detector ID = " << subdet << ": file " << __FILE__ << " line " << __LINE__ ;
              throw GaudiException( errorMsg.str(), "CellID not matched", StatusCode::FAILURE );  
            }

            // add tolerence 
            strip_length_mm = strip_length_mm * (1.0 + _striplength_tolerance);
	    try{
	      edm4hep::TrackerHit spacePoint = createSpacePoint( hitFront, hitBack, strip_length_mm);

	      //UTIL::CellIDEncoder<TrackerHitImpl> cellid_encoder( UTIL::ILDCellID0::encoder_string , spCol );
              //cellid_encoder.setValue( cellID0 ); //give the new hit, the CellID0 of the front hit
              //cellid_encoder.setCellID( spacePoint ) ;
	      spacePoint.setCellID(cellID0);

              // store the hits it's composed of:
              spacePoint.addToRawHits( hitFront.getObjectID() );
              spacePoint.addToRawHits( hitBack.getObjectID() );

              spacePoint.setType( UTIL::set_bit( spacePoint.getType() ,  UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ) ) ;

              spCol->push_back( spacePoint );
	      debug() << "Hit accepted id = " << spacePoint.id() << " cellID = " << spacePoint.getCellID() << endmsg;
	      //debug() << "push_back space point's id=" << spCol->at(spCol->size()-1).id() << endmsg;
              createdSpacePoints++;
              
              ///////////////////////////////
              // make the relations
              if( simHitsFront.size() == 1 ){
		edm4hep::ConstSimTrackerHit& simHit = simHitsFront[0];
		edm4hep::MCRecoTrackerAssociation spAss = relCol->create();
		spAss.setRec(spacePoint);
		spAss.setSim(simHit);
		spAss.setWeight( 0.5 );
              }
              if( simHitsBack.size() == 1 ){
		edm4hep::ConstSimTrackerHit& simHit = simHitsBack[0];
		edm4hep::MCRecoTrackerAssociation spAss = relCol->create();
                spAss.setRec(spacePoint);
		spAss.setSim(simHit);
                spAss.setWeight( 0.5 );
	      }
	    }
	    catch(std::runtime_error& e){
	      if ( ghost_hit == true ) {
                debug() << "Ghosthit correctly rejected" << endmsg;
              }
	      else {
                debug() << "True hit rejected!" << endmsg;
              }
            }
          }
        }
      }
    }
    
    //if(spCol->size()!=0) _outSPColHdl.put(spCol);
    //else                 delete spCol;
    //if(relCol->size()!=0) _outSPAssColHdl.put(relCol);
    //else                  delete relCol;
    
    debug() << "\nCreated " << createdSpacePoints << " space points ( raw strip hits: " << rawStripHits << ")\n"
	    << "  There were " << rawStripHits << " strip hits available, giving " 
	    << possibleSpacePoints << " possible space points" << endmsg;
    debug() << "  " << _nStripsTooParallel << " space points couldn't be created, because the strips were too parallel\n"
	    << "  " << _nPlanesNotParallel << " space points couldn't be created, because the planes of the measurement surfaces where not parallel enough\n"
	    << "  " << _nOutOfBoundary     << " space points couldn't be created, because the result was outside the sensor boundary\n" << endmsg; 
  }

  _nEvt ++ ;
  return sc;
}

StatusCode SpacePointBuilderAlg::finalize(){
  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

edm4hep::TrackerHit SpacePointBuilderAlg::createSpacePoint( edm4hep::ConstTrackerHit a , edm4hep::ConstTrackerHit b, double stripLength ){
  
  const edm4hep::Vector3d& pa = a.getPosition();
  double xa = pa[0];
  double ya = pa[1];
  double za = pa[2];
  CLHEP::Hep3Vector PA( xa,ya,za );
  //double du_a = a->getdU();  
  float du_a = a.getCovMatrix(2);
  
  gear::MeasurementSurface const* msA = _GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( a.getCellID() );
  gear::CartesianCoordinateSystem* ccsA = dynamic_cast< gear::CartesianCoordinateSystem* >( msA->getCoordinateSystem() );
 
  CLHEP::Hep3Vector WA = ccsA->getLocalZAxis(); // the vector W of the local coordinate system the measurement surface has
  CLHEP::Hep3Vector VA = ccsA->getLocalYAxis(); // the vector V of the local coordinate system the measurement surface has
  CLHEP::Hep3Vector UA = ccsA->getLocalXAxis(); // the vector U of the local coordinate system the measurement surface has
  
  const edm4hep::Vector3d& pb = b.getPosition();
  double xb = pb[0];
  double yb = pb[1];
  double zb = pb[2];
  CLHEP::Hep3Vector PB( xb,yb,zb );
  //double du_b = b->getdU();  
  float du_b = b.getCovMatrix(2);

  gear::MeasurementSurface const* msB = _GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( b.getCellID() );
  gear::CartesianCoordinateSystem* ccsB = dynamic_cast< gear::CartesianCoordinateSystem* >( msB->getCoordinateSystem() );
  CLHEP::Hep3Vector WB = ccsB->getLocalZAxis(); // the vector W of the local coordinate system the measurement surface has
  CLHEP::Hep3Vector VB = ccsB->getLocalYAxis(); // the vector V of the local coordinate system the measurement surface has
  CLHEP::Hep3Vector UB = ccsB->getLocalXAxis(); // the vector U of the local coordinate system the measurement surface has
  
  debug() << "\t ( " << xa << " " << ya << " " << za << " ) <--> ( " << xb << " " << yb << " " << zb << " )" << endmsg;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // First: check if the two measurement surfaces are parallel (i.e. the w are parallel or antiparallel)
  double angle = fabs(WB.angle(WA));
  double angleMax = 1.*M_PI/180.;
  if(( angle > angleMax )&&( angle < M_PI-angleMax )){
    _nPlanesNotParallel++;
    debug() << "\tThe planes of the measurement surfaces are not parallel enough, the angle between the W vectors is " << angle
	    << " where the angle has to be smaller than " << angleMax << " or bigger than " << M_PI-angleMax << endmsg;
    throw std::runtime_error(""); //calculate the xing point and if that fails don't create a spacepoint
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Next: check if the angle between the strips is not 0
  angle = fabs(VB.angle(VA));
  double angleMin= 1.*M_PI/180.;
  if(( angle < angleMin )||( angle > M_PI-angleMin )){
    _nStripsTooParallel++;
    debug() << "\tThe strips (V vectors) of the measurement surfaces are too parallel, the angle between the V vectors is " << angle
	    << " where the angle has to be between " << angleMax << " or bigger than " << M_PI-angleMin << endmsg;
    throw std::runtime_error(""); //calculate the xing point and if that fails don't create a spacepoint
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Next we want to calculate the crossing point.
  
  CLHEP::Hep3Vector point;

//  calculatePointBetweenTwoLines( PA, VA, PB, VB, point );
//  
//  // we want to set the space point on the surface of the hit closest to the IP
//  if (PA.mag2() < PB.mag2()) {
//    calculatePointBetweenTwoLines( PA, VA, PB, VB, point );
//  } else {
//    calculatePointBetweenTwoLines( PB, VB, PA, VA, point );
//  }
//  
//
//  
//  streamlog_out( DEBUG2 ) << "\tStandard: Position of space point (global) : ( " << point.x() << " " << point.y() << " " << point.z() << " )\n";
   
  CLHEP::Hep3Vector vertex(0.,0.,0.);
  CLHEP::Hep3Vector L1 = ccsA->getLocalPoint(PA);
  CLHEP::Hep3Vector L2 = ccsB->getLocalPoint(PB);

  debug() << " L1 = " << L1 << endmsg;
  debug() << " L2 = " << L2 << endmsg;
      
  L1.setY(-stripLength/2.0);
  CLHEP::Hep3Vector SL1 = L1;
  L1.setY( stripLength/2.0);
  CLHEP::Hep3Vector EL1 = L1;
  
  
  L2.setY(-stripLength/2.0);
  CLHEP::Hep3Vector SL2 = L2;
  L2.setY( stripLength/2.0);
  CLHEP::Hep3Vector EL2 = L2;
  
  
  CLHEP::Hep3Vector S1 = ccsA->getGlobalPoint(SL1);
  CLHEP::Hep3Vector E1 = ccsA->getGlobalPoint(EL1);

  CLHEP::Hep3Vector S2 = ccsB->getGlobalPoint(SL2);
  CLHEP::Hep3Vector E2 = ccsB->getGlobalPoint(EL2);

  debug() << " stripLength = " << stripLength << endmsg;
  debug() << " S1 = " << S1 << endmsg;
  debug() << " E1 = " << E1 << endmsg;
  debug() << " S2 = " << S2 << endmsg;
  debug() << " E2 = " << E2 << endmsg;

  point.set(0.0, 0.0, 0.0);
    
  int valid_intersection = calculatePointBetweenTwoLines_UsingVertex( S1, E1, S2, E2, vertex, point );
  
  if (valid_intersection != 0) {
    debug() << "\tNo valid intersection for lines" << endmsg;
    throw std::runtime_error("");
  }
  
  debug() << "\tVertex: Position of space point (global) : ( " << point.x() << " " << point.y() << " " << point.z() << " )" << endmsg;
  
  // Check if the new hit is within the boundaries
  CLHEP::Hep3Vector localPointA = ccsA->getLocalPoint(point);
  localPointA.setZ(0.);// we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  
  CLHEP::Hep3Vector localPointB = ccsB->getLocalPoint(point);
  localPointB.setZ(0.);// we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  
  if( !msA->isLocalInBoundary( localPointA ) ){
    _nOutOfBoundary++;
    debug() << "\tSpacePoint position lies outside the boundary of the first layer: local coordinates are ( " 
	    << localPointA.x() << " " << localPointA.y() << " " << localPointA.z() << " )" << endmsg;
    throw std::runtime_error("");
  }
  if( !msB->isLocalInBoundary( localPointB ) ){
    _nOutOfBoundary++;
    debug() << "\tSecond hit is out of boundary: local coordinates are ( " 
	    << localPointB.x() << " " << localPointB.y() << " " << localPointB.z() << " )" << endmsg;
    throw std::runtime_error("");
  }
  
  //Create the new TrackerHit
  edm4hep::TrackerHit spacePoint;// = new edm4hep::TrackerHit();
  
  edm4hep::Vector3d pos(point.x(), point.y(), point.z());
  spacePoint.setPosition(pos) ;
  
  // set error treating the strips as stereo with equal and opposite rotation -- for reference see Karimaki NIM A 374 p367-370

  // first calculate the covariance matrix in the cartisian coordinate system defined by the sensor 
  // here we assume that du is the same for both sides
  
  if( fabs(du_a - du_b) > 1.0e-06 ){
    error() << "\tThe measurement errors of the two 1D hits must be equal " << endmsg;    
    assert( fabs(du_a - du_b) > 1.0e-06 == false );
    throw std::runtime_error(""); //measurement errors are not equal don't create a spacepoint
  }
   
  double du2 = du_a*du_a;
  
  // rotate the strip system back to double-layer wafer system
  CLHEP::Hep3Vector u_sensor = UA + UB;
  CLHEP::Hep3Vector v_sensor = VA + VB;
  CLHEP::Hep3Vector w_sensor = WA + WB;
  
  CLHEP::HepRotation rot_sensor( u_sensor, v_sensor, w_sensor );
  CLHEP::HepMatrix rot_sensor_matrix;
  rot_sensor_matrix = rot_sensor;
  
  double cos2_alpha = VA.cos2Theta(v_sensor) ; // alpha = strip angle   
  double sin2_alpha = 1 - cos2_alpha ; 
  
  CLHEP::HepSymMatrix cov_plane(3,0); // u,v,w
  
  cov_plane(1,1) = (0.5 * du2) / cos2_alpha;
  cov_plane(2,2) = (0.5 * du2) / sin2_alpha;
  
  debug() << "cov_plane  = " << cov_plane << endmsg;  
  debug() << "strip_angle = " << VA.angle(VB)/(M_PI/180) / 2.0 << " degrees " << endmsg;
  
  CLHEP::HepSymMatrix cov_xyz= cov_plane.similarity(rot_sensor_matrix);
  
  debug() << "cov_xyz  = " << cov_xyz << endmsg;
  
  std::array<float, 6> cov;
  //EVENT::FloatVec cov( 9 )  ; 
  int icov = 0 ;
  for(int irow=0; irow<3; ++irow ){
    for(int jcol=0; jcol<irow+1; ++jcol){
      //      std::cout << "row = " << irow << " col = " << jcol << std::endl ;
      cov[icov] = cov_xyz[irow][jcol] ;
//      std::cout << "cov["<< icov << "] = " << cov[icov] << std::endl ;
      ++icov ;
    }
  }
  
  spacePoint.setCovMatrix(cov);
    
  return spacePoint;
}
/*
TrackerHitImpl* SpacePointBuilderAlg::createSpacePointOld( TrackerHitPlane* a , TrackerHitPlane* b ){
  
  streamlog_out( DEBUG2 ) << "\t OLD OLD OLD OLD\n";
  
  
  const double* p1 = a.getPosition();
  double x1 = p1[0];
  double y1 = p1[1];
  double z1 = p1[2];
  const float* v1 = a->getV();
  float ex1 = cos( v1[1] ) * sin( v1[0] ); 
  float ey1 = sin( v1[1] ) * sin( v1[0] );
  
  const double* p2 = b.getPosition();
  double x2 = p2[0];
  double y2 = p2[1];
  double z2 = p2[2];
  const float* v2 = b->getV();
  float ex2 = cos( v2[1] ) * sin( v2[0] ); 
  float ey2 = sin( v2[1] ) * sin( v2[0] );
  
  streamlog_out( DEBUG2 ) << "\t ( " << x1 << " " << y1 << " " << z1 << " ) <--> ( " << x2 << " " << y2 << " " << z2 << " )\n";
  
  double x=0.;
  double y=0.;
  
  if ( calculateXingPoint( x1, y1, ex1, ey1, x2, y2, ex2, ey2, x, y ) != 0 ){
    
    _nStripsTooParallel++;
    streamlog_out( DEBUG2 ) << "\tStrips too parallel\n\n";
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
  
  }
  
  double z= (z1 + z2)/2.;
  
  streamlog_out( DEBUG2 ) << "\tPosition of space point (global) : ( " << x << " " << y << " " << z << " )\n";
  
  // Check if the new hit is within the boundary
  CLHEP::Hep3Vector globalPoint(x,y,z);
  //  gear::MeasurementSurface* ms = gear::MeasurementSurfaceStore::Instance().GetMeasurementSurface( a.getCellID() );
  CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
  localPoint.setZ( 0. ); // we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  if( !ms->isLocalInBoundary( localPoint ) ){
    
    _nOutOfBoundary++;
    streamlog_out( DEBUG2 ) << "\tHit is out of boundary: local coordinates are ( " 
      << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )\n\n";
    
    return NULL;
    
  }
  
  
  //Create the new TrackerHit
  TrackerHitImpl* spacePoint = new TrackerHitImpl();
  
  double pos[3] = {x,y,z};
  spacePoint->setPosition(  pos  ) ;
  
  streamlog_out( DEBUG2 ) << "\tHit accepted\n\n";

  return spacePoint;
  
}
*/



int SpacePointBuilderAlg::calculatePointBetweenTwoLines_UsingVertex( 
                                                const CLHEP::Hep3Vector& PA, 
                                                const CLHEP::Hep3Vector& PB, 
                                                const CLHEP::Hep3Vector& PC, 
                                                const CLHEP::Hep3Vector& PD,
                                                const CLHEP::Hep3Vector& Vertex,
                                                CLHEP::Hep3Vector& point){

  
  // A general point on the line joining point PA to point PB is
  // x, where 2*x=(1+m)*PA + (1-m)*PB. Similarly for 2*y=(1+n)*PC + (1-n)*PD.
  // Suppose that v is the vertex. Requiring that the two 'general
  // points' lie on a straight through v means that the vector x-v is a 
  // multiple of y-v. This condition fixes the parameters m and n.
  // We then return the 'space-point' x, supposed to be the layer containing PA and PB. 
  // We require that -1<m<1, otherwise x lies 
  // outside the segment PA to PB; and similarly for n.
  
  bool ok = true;
  
//  streamlog_out( DEBUG1 ) << " Vertex = " << Vertex << std::endl; 
//  
//  streamlog_out( DEBUG1 ) << " PA = " << PA << std::endl;
//  streamlog_out( DEBUG1 ) << " PB = " << PB << std::endl;
//  streamlog_out( DEBUG1 ) << " PC = " << PC << std::endl;
//  streamlog_out( DEBUG1 ) << " PD = " << PD << std::endl;
  
  CLHEP::Hep3Vector VAB(PA-PB);
  CLHEP::Hep3Vector VCD(PC-PD);

//  streamlog_out( DEBUG1 ) << " VAB = " << VAB << std::endl;
//  streamlog_out( DEBUG1 ) << " VCD = " << VCD << std::endl;
  
  CLHEP::Hep3Vector  s(PA+PB-2*Vertex);   // twice the vector from vertex to midpoint
  CLHEP::Hep3Vector  t(PC+PD-2*Vertex);   // twice the vector from vertex to midpoint

  CLHEP::Hep3Vector  qs(VAB.cross(s));  
  CLHEP::Hep3Vector  rt(VCD.cross(t));  

//  streamlog_out( DEBUG1 ) << " s = " << s << std::endl;
//  streamlog_out( DEBUG1 ) << " t = " << t << std::endl;
//  streamlog_out( DEBUG1 ) << " qs = " << qs << std::endl;
//  streamlog_out( DEBUG1 ) << " rt = " << rt << std::endl;
  
  
  double m = (-(s*rt)/(VAB*rt)); // ratio for first line
    
  double limit = 1.0;
  
  if (m>limit || m<-1.*limit) {
    debug() << "m' = " << m << endmsg;
    ok = false;
  }
  else {
    double n = (-(t*qs)/(VCD*qs)); // ratio for second line
    if (n>limit || n<-1.*limit) {
      debug() << "n' = " << n << endmsg;
      ok = false;
    }
  }
  
  if (ok) {
    point = 0.5*(PA + PB + m*VAB);
  }
  
  return ok ? 0 : 1;
}



int SpacePointBuilderAlg::calculatePointBetweenTwoLines( const CLHEP::Hep3Vector& P1, const CLHEP::Hep3Vector& V1, const CLHEP::Hep3Vector& P2, const CLHEP::Hep3Vector& V2, CLHEP::Hep3Vector& point ){
  
  // Richgungsvektor normal auf die anderen beiden:
  CLHEP::Hep3Vector n = V1.cross( V2 );
  
  // Now we want to rotate into a coordinate system, where n is parallel to the z axis
  // For this: first set phi to 0
  // then: set theta to 0 (we set phi to 0 first, so we can then rotate arount the y axis)
  CLHEP::HepRotation rot;
  rot.rotateZ( -n.phi() );
  CLHEP::Hep3Vector nPrime = rot * n; //now the phi of nPrime should be 0
  debug() << "phi of n'' = " << nPrime.phi() << " (it should be 0!!!)\n";
  rot.rotateY( -n.theta() );
  nPrime = rot * n;
  debug() << "phi of n' = " << nPrime.phi() << " (it should be 0!!!)\n";
  debug() << "theta of n' = " << nPrime.theta() <<  " (it should be 0!!!)" << endmsg;
  
  // Now rotate all the vectors and points into this coordinatesystem.
  CLHEP::Hep3Vector P1prime = rot * P1;
  CLHEP::Hep3Vector V1prime = rot * V1;
  CLHEP::Hep3Vector P2prime = rot * P2;
  CLHEP::Hep3Vector V2prime = rot * V2;
  
  // What is the gain of rotating into this system?
  // A: 
  double x;
  double y;
  int res = calculateXingPoint( P1prime.x(), P1prime.y(), V1prime.x(), V1prime.y(), P2prime.x(), P2prime.y(), V2prime.x(), V2prime.y(), x, y );
  
  if ( res != 0 ) return 1;
  
  point.setX( x );
  point.setY( y );
  point.setZ( (P1prime.z() + P2prime.z())/2. );
  
  // Now transform back to the global coordinates
  point = rot.inverse() * point;
  
  return 0;
}


int SpacePointBuilderAlg::calculateXingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y ){
  float a = (x1*ey1 - y1*ex1) - (x2*ey1 - y2*ex1);
  float b = ex2*ey1 - ex1*ey2;

  const float epsilon = 0.00001;

  if( fabs(b) < epsilon ) return 1; // if b==0 the two directions e1 and e2 are parallel and there is no crossing!

  float t = a/b;

  x = x2 + t*ex2;
  y = y2 + t*ey2;

  return 0;
}
 
 
std::vector< int > SpacePointBuilderAlg::getCellID0sAtBack( int cellID0 ){
  std::vector< int > back;
  
  UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
  int subdet = cellID[ UTIL::ILDCellID0::subdet ] ;
  
  if( subdet == UTIL::ILDDetID::FTD ) return getCellID0sAtBackOfFTD( cellID0 );
  if( subdet == UTIL::ILDDetID::SIT ) return getCellID0sAtBackOfSIT( cellID0 );
  if( subdet == UTIL::ILDDetID::SET ) return getCellID0sAtBackOfSET( cellID0 );
  
  return back;
}


std::vector< int > SpacePointBuilderAlg::getCellID0sAtBackOfFTD( int cellID0 ){
  std::vector< int > back;
  
  const gear::FTDLayerLayout& ftdLayers = _GEAR->getFTDParameters().getFTDLayerLayout();
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
//   int side   = cellID[ UTIL::ILDCellID0::side ];
//   int module = cellID[ UTIL::ILDCellID0::module ];
  int sensor = cellID[ UTIL::ILDCellID0::sensor ];
  int layer  = cellID[ UTIL::ILDCellID0::layer ];
  
  //check if sensor is in front
  if(( ftdLayers.isDoubleSided( layer ) ) && ( sensor <= ftdLayers.getNSensors( layer ) / 2 ) ){
    cellID[ UTIL::ILDCellID0::sensor ] = sensor + ftdLayers.getNSensors( layer ) / 2; 
    // it is assumed (according to current gear and mokka), that sensors 1 until n/2 will be on front
    // and sensor n/2 + 1 until n are at the back
    // so the sensor x, will have sensor x+n/2 at the back
    back.push_back( cellID.lowWord() );
  }
  return back;
}

std::vector< int > SpacePointBuilderAlg::getCellID0sAtBackOfSIT( int cellID0 ){
  std::vector< int > back;
  
//   const gear::ZPlanarLayerLayout& sitLayout = _GEAR->getSITParameters().getZPlanarLayerLayout();
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
//   int side   = cellID[ UTIL::ILDCellID0::side ];
//   int module = cellID[ UTIL::ILDCellID0::module ];
//   int sensor = cellID[ UTIL::ILDCellID0::sensor ];
  int layer  = cellID[ UTIL::ILDCellID0::layer ];
  
  //check if sensor is in front
  if( layer%2 == 0 ){ // even layers are front sensors
    cellID[ UTIL::ILDCellID0::layer ] = layer + 1; 
    // it is assumed that the even layers are the front layers
    // and the following odd ones the back layers
    back.push_back( cellID.lowWord() );
  }
  return back;
}

std::vector< int > SpacePointBuilderAlg::getCellID0sAtBackOfSET( int cellID0 ){
  std::vector< int > back;
  
//   const gear::ZPlanarLayerLayout& setLayout = _GEAR->getSETParameters().getZPlanarLayerLayout();
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
  //   int side   = cellID[ UTIL::ILDCellID0::side ];
  //   int module = cellID[ UTIL::ILDCellID0::module ];
  //   int sensor = cellID[ UTIL::ILDCellID0::sensor ];
  int layer  = cellID[ UTIL::ILDCellID0::layer ];
  
  //check if sensor is in front
  if( layer%2 == 0 ){ // even layers are front sensors
    cellID[ UTIL::ILDCellID0::layer ] = layer + 1; 
    // it is assumed that the even layers are the front layers
    // and the following odd ones the back layers
    back.push_back( cellID.lowWord() );
  }
  return back;
}

std::string SpacePointBuilderAlg::getCellID0Info( int cellID0 ){
  std::stringstream s;
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );

  int subdet = cellID[ UTIL::ILDCellID0::subdet ] ;
  int side   = cellID[ UTIL::ILDCellID0::side ];
  int module = cellID[ UTIL::ILDCellID0::module ];
  int sensor = cellID[ UTIL::ILDCellID0::sensor ];
  int layer  = cellID[ UTIL::ILDCellID0::layer ];
  
  s << "(su" << subdet << ",si" << side << ",la" << layer << ",mo" << module << ",se" << sensor << ")";
  
  return s.str();
}



