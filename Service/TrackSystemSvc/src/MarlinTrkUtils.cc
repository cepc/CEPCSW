#include "TrackSystemSvc/MarlinTrkUtils.h"

#include <vector>
#include <algorithm>

#include "TrackSystemSvc/IMarlinTrack.h"
#include "TrackSystemSvc/HelixTrack.h"

#include "DataHelper/Navigation.h"
//#include "lcio.h"
//#include <IMPL/TrackImpl.h>
//#include <IMPL/TrackStateImpl.h>
//#include <EVENT/TrackerHit.h>
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/Track.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

//#include "streamlog/streamlog.h"

#include "TMatrixD.h"

#define MIN_NDF 6

namespace MarlinTrk {
  

//  // Check if a square matrix is Positive Definite 
//  bool Matrix_Is_Positive_Definite(const EVENT::FloatVec& matrix){
//    
//    std::cout << "\n MarlinTrk::Matrix_Is_Positive_Definite(EVENT::FloatVec& matrix): " << std::endl;
//    
//    int icol,irow;
//    
//    int nrows = 5;
//    
//    TMatrixD cov(nrows,nrows) ; 
//    
//    bool matrix_is_positive_definite = true;
//    
//    int icov = 0;
//    for(int irow=0; irow<nrows; ++irow ){
//      for(int jcol=0; jcol<irow+1; ++jcol){
//        cov(irow,jcol) = matrix[icov];
//        cov(jcol,irow) = matrix[icov];
////        std::cout << " " << matrix[icov] ;
//        ++icov ;
//      }
////      std::cout << std::endl;
//    }
//    
////    cov.Print();
//    
//    double *pU = cov.GetMatrixArray();
//    
//    for (icol = 0; icol < nrows; icol++) {
//      const int rowOff = icol * nrows;
//      
//      //Compute fU(j,j) and test for non-positive-definiteness.
//      double ujj = pU[rowOff+icol];
//      double diagonal = ujj;
////      std::cout << "ERROR: diagonal = " << diagonal << std::endl;
//
//      for (irow = 0; irow < icol; irow++) {
//        const int pos_ij = irow*nrows+icol;
//        std::cout << " " << pU[pos_ij] ;
//        ujj -= pU[pos_ij]*pU[pos_ij];
//      }
//      std::cout  << " " << diagonal << std::endl;
//
//      
//      if (ujj <= 0) {
//        matrix_is_positive_definite = false;
//      }
//    }
//    
//      std::cout  << std::endl;
//    
//    if ( matrix_is_positive_definite == false ) {
//      std::cout << "******************************************************" << std::endl;
//      std::cout << "** ERROR:  matrix shown not to be positive definite **" << std::endl;
//      std::cout << "******************************************************" << std::endl;
//    }
//    
//    return matrix_is_positive_definite;
//    
//  }
  
  
  
  int createTrackStateAtCaloFace( IMarlinTrack* marlinTrk, edm4hep::TrackState* track, edm4hep::ConstTrackerHit trkhit, bool tanL_is_positive );
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<edm4hep::ConstTrackerHit>& hit_list, edm4hep::Track* track, bool fit_backwards, const std::array<float,15>& initial_cov_for_prefit, float bfield_z, double maxChi2Increment){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs ;
    
    if( track == 0 ){
      throw std::runtime_error( "MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL " ) ;
    }
    
    int return_error = 0;
    ///////////////////////////////////////////////////////
    // produce prefit parameters
    ///////////////////////////////////////////////////////
    
    edm4hep::TrackState pre_fit ;
    
    //std::cout << "debug:=====================before createPrefit" << std::endl;
    return_error = createPrefit(hit_list, &pre_fit, bfield_z, fit_backwards);
    //std::cout << "debug:=====================after createPrefit return code=" << return_error << std::endl;
    pre_fit.covMatrix = initial_cov_for_prefit;

    ///////////////////////////////////////////////////////
    // use prefit parameters to produce Finalised track
    ///////////////////////////////////////////////////////
    
    if( return_error == 0 ) {
      
      return_error = createFinalisedLCIOTrack( marlinTrk, hit_list, track, fit_backwards, &pre_fit, bfield_z, maxChi2Increment);
      
    } else {
      //std::cout << "MarlinTrk::createFinalisedLCIOTrack : Prefit failed error = " << return_error << std::endl;
    }
    return return_error;
  }
  
  int createFinalisedLCIOTrack( IMarlinTrack* marlinTrk, std::vector<edm4hep::ConstTrackerHit>& hit_list, edm4hep::Track* track, bool fit_backwards, edm4hep::TrackState* pre_fit, float bfield_z, double maxChi2Increment){
    
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs ;
    
    if( track == 0 ){
      throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ");
    }
    
    if( pre_fit == 0 ){
      throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ");
    }
    
    
    int fit_status = createFit(hit_list, marlinTrk, pre_fit, bfield_z, fit_backwards, maxChi2Increment);
    
    if( fit_status != IMarlinTrack::success ){ 
      
      //std::cout << "MarlinTrk::createFinalisedLCIOTrack fit failed: fit_status = " << fit_status << std::endl; 
      
      return fit_status;
    } 
    
    int error = finaliseLCIOTrack(marlinTrk, track, hit_list);
        
    return error;
  }
  
  
  
  
  int createFit( std::vector<edm4hep::ConstTrackerHit>& hit_list, IMarlinTrack* marlinTrk, edm4hep::TrackState* pre_fit, float bfield_z, bool fit_backwards, double maxChi2Increment){
    
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs;
    
    if( marlinTrk == 0 ){
      throw std::runtime_error("MarlinTrk::createFit: IMarlinTrack == NULL ");
    }
    
    if( pre_fit == 0 ){
      throw std::runtime_error("MarlinTrk::createFit: TrackStateImpl == NULL ");
    }
    
    int return_error = 0;
    
    
//    ///////////////////////////////////////////////////////
//    // check that the prefit has the reference point at the correct location 
//    ///////////////////////////////////////////////////////
//    
//    if (( fit_backwards == IMarlinTrack::backward && pre_fit->getLocation() != lcio::TrackState::AtLastHit ) 
//        ||  
//        ( fit_backwards == IMarlinTrack::forward && pre_fit->getLocation() != lcio::TrackState::AtFirstHit )) {            
//      std::stringstream ss ;
//      
//      ss << "MarlinTrk::createFinalisedLCIOTrack track state must be set at either first or last hit. Location = ";
//      ss << pre_fit->getLocation();
//      
//      throw EVENT::Exception( ss.str() );
//      
//    } 
    
    ///////////////////////////////////////////////////////
    // add hits to IMarlinTrk  
    ///////////////////////////////////////////////////////
    
    std::vector<edm4hep::ConstTrackerHit>::iterator it = hit_list.begin();
    
    //  start by trying to add the hits to the track we want to finally use. 
    //std::cout << "MarlinTrk::createFit Start Fit: AddHits: number of hits to fit " << hit_list.size() << std::endl;
    
    std::vector<edm4hep::ConstTrackerHit> added_hits;
    unsigned int ndof_added = 0;
    
    for( it = hit_list.begin() ; it != hit_list.end() ; ++it ) {
      
      edm4hep::ConstTrackerHit trkHit = *it;
      bool isSuccessful = false; 
      //std::cout << "debug: TrackerHit pointer " << trkHit << std::endl;
      if( UTIL::BitSet32( trkHit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint        
        //Split it up and add both hits to the MarlinTrk
        //const EVENT::LCObjectVec rawObjects = trkHit->getRawHits();
	//std::cout << "space point is not still valid! pelease wait updating..." <<std::endl;
	//exit(1);
        int nRawHit = trkHit.rawHits_size();
        for( unsigned k=0; k< nRawHit; k++ ){
          edm4hep::ConstTrackerHit rawHit = Navigation::Instance()->GetTrackerHit(trkHit.getRawHits(k));
	  if( marlinTrk->addHit( rawHit ) == IMarlinTrack::success ){
	    isSuccessful = true; //if at least one hit from the spacepoint gets added
            ++ndof_added;
//            streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;
          }
        }
      }
      else { // normal non composite hit
        if (marlinTrk->addHit( trkHit ) == IMarlinTrack::success ) {
          isSuccessful = true;
          ndof_added += 2;
//          streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;          
        }
      }
      
      if (isSuccessful) {
        added_hits.push_back(trkHit);
      }        
      else{
	//std::cout << "MarlinTrkUtils::createFit Hit " << it - hit_list.begin() << " Dropped " << std::endl;
      }
      
    }
      
    if( ndof_added < MIN_NDF ) {
      //streamlog_out(DEBUG2) << "MarlinTrk::createFit : Cannot fit less with less than " << MIN_NDF << " degrees of freedom. Number of hits =  " << added_hits.size() << " ndof = " << ndof_added << std::endl;
      return IMarlinTrack::bad_intputs;
    }
      
    
    
    ///////////////////////////////////////////////////////
    // set the initial track parameters  
    ///////////////////////////////////////////////////////
    
    return_error = marlinTrk->initialise( *pre_fit, bfield_z, IMarlinTrack::backward ) ;
    if (return_error != IMarlinTrack::success) {
      
      //streamlog_out(DEBUG5) << "MarlinTrk::createFit Initialisation of track fit failed with error : " << return_error << std::endl;
      
      return return_error;
      
    }
    
    ///////////////////////////////////////////////////////
    // try fit and return error
    ///////////////////////////////////////////////////////
    int status = marlinTrk->fit(maxChi2Increment);
    //std::cout << "debug:===================createFit " << status << std::endl;
    
    return status;
    
  }
  
  
  
  int createPrefit( std::vector<edm4hep::ConstTrackerHit>& hit_list, edm4hep::TrackState* pre_fit, float bfield_z, bool fit_backwards){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if ( hit_list.empty() ) return IMarlinTrack::bad_intputs ;
    
    if( pre_fit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackStateImpl == NULL ")  ) ;
    }
    
    ///////////////////////////////////////////////////////
    // loop over all the hits and create a list consisting only 2D hits 
    ///////////////////////////////////////////////////////
    
    std::vector<edm4hep::ConstTrackerHit> twoD_hits;
    
    for (unsigned ihit=0; ihit < hit_list.size(); ++ihit) {
      
      // check if this a space point or 2D hit 
      if(UTIL::BitSet32( hit_list[ihit].getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] == false ){
        // then add to the list 
        twoD_hits.push_back(hit_list[ihit]);
        
      }
    }
    
    ///////////////////////////////////////////////////////
    // check that there are enough 2-D hits to create a helix 
    ///////////////////////////////////////////////////////
    
    if (twoD_hits.size() < 3) { // no chance to initialise print warning and return
      //streamlog_out(WARNING) << "MarlinTrk::createFinalisedLCIOTrack Cannot create helix from less than 3 2-D hits" << std::endl;
      return IMarlinTrack::bad_intputs;
    }
    
    ///////////////////////////////////////////////////////
    // make a helix from 3 hits to get a trackstate
    ///////////////////////////////////////////////////////
    
    // SJA:FIXME: this may not be the optimal 3 hits to take in certain cases where the 3 hits are not well spread over the track length 
    const edm4hep::Vector3d& x1 = twoD_hits[0].getPosition();
    const edm4hep::Vector3d& x2 = twoD_hits[ twoD_hits.size()/2 ].getPosition();
    const edm4hep::Vector3d& x3 = twoD_hits.back().getPosition();
    
    HelixTrack helixTrack( x1, x2, x3, bfield_z, HelixTrack::forwards );

    if ( fit_backwards == IMarlinTrack::backward ) {
      pre_fit->location = MarlinTrk::Location::AtLastHit;
      helixTrack.moveRefPoint(hit_list.back().getPosition()[0], hit_list.back().getPosition()[1], hit_list.back().getPosition()[2]);      
    } else {
      pre_fit->location = MarlinTrk::Location::AtFirstHit;
      helixTrack.moveRefPoint(hit_list.front().getPosition()[0], hit_list.front().getPosition()[1], hit_list.front().getPosition()[2]);      
    }
    
    
    const float referencePoint[3] = { helixTrack.getRefPointX() , helixTrack.getRefPointY() , helixTrack.getRefPointZ() };
    
    pre_fit->D0 = helixTrack.getD0();
    pre_fit->phi = helixTrack.getPhi0();
    pre_fit->omega = helixTrack.getOmega();
    pre_fit->Z0 = helixTrack.getZ0();
    pre_fit->tanLambda = helixTrack.getTanLambda();
    
    pre_fit->referencePoint = referencePoint;
    
    return IMarlinTrack::success;
    
  }
  
  int finaliseLCIOTrack( IMarlinTrack* marlintrk, edm4hep::Track* track, std::vector<edm4hep::ConstTrackerHit>& hit_list, edm4hep::TrackState* atLastHit, edm4hep::TrackState* atCaloFace){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( marlintrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: IMarlinTrack == NULL ")  ) ;
    }
    
    if( track == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: TrackImpl == NULL ")  ) ;
    }
    
    if( atCaloFace && atLastHit == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: atLastHit == NULL ")  ) ;
    }

    if( atLastHit && atCaloFace == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::finaliseLCIOTrack: atCaloFace == NULL ")  ) ;
    }

    
    
    ///////////////////////////////////////////////////////
    // error to return if any
    ///////////////////////////////////////////////////////
    int return_error = 0;
    
    int ndf = 0;
    double chi2 = -DBL_MAX;
    
    /////////////////////////////////////////////////////////////
    // First check NDF to see if it make any sense to continue.
    // The track will be dropped if the NDF is less than 0 
    /////////////////////////////////////////////////////////////
    
    return_error = marlintrk->getNDF(ndf);

    if ( return_error != IMarlinTrack::success) {
      //streamlog_out(DEBUG3) << "MarlinTrk::finaliseLCIOTrack: getNDF returns " << return_error << std::endl;
      return return_error;
    } else if( ndf < 0 ) {
      //streamlog_out(DEBUG2) << "MarlinTrk::finaliseLCIOTrack: number of degrees of freedom less than 0 track dropped : NDF = " << ndf << std::endl;
      return IMarlinTrack::error;
    } else {
      //streamlog_out(DEBUG1) << "MarlinTrk::finaliseLCIOTrack: NDF = " << ndf << std::endl;
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // get the list of hits used in the fit
    // add these to the track, add spacepoints as long as at least on strip hit is used.  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hits_in_fit;
    std::vector<std::pair<edm4hep::ConstTrackerHit, double> > outliers;
    std::vector<edm4hep::ConstTrackerHit> used_hits;
        
    hits_in_fit.reserve(300);
    outliers.reserve(300);
    
    marlintrk->getHitsInFit(hits_in_fit);
    marlintrk->getOutliers(outliers);
    
    ///////////////////////////////////////////////
    // now loop over the hits provided for fitting 
    // we do this so that the hits are added in the
    // order in which they have been fitted
    ///////////////////////////////////////////////
    
    for ( unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {
      
      edm4hep::ConstTrackerHit trkHit = hit_list[ihit];
      
      if( UTIL::BitSet32( trkHit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
	//std::cout << "Error: space point is not still valid! pelease wait updating..." <<std::endl;
        //exit(1);
	// get strip hits 
        int nRawHit = trkHit.rawHits_size();
        for( unsigned k=0; k< nRawHit; k++ ){
	  edm4hep::ConstTrackerHit rawHit = Navigation::Instance()->GetTrackerHit(trkHit.getRawHits(k));
	  bool is_outlier = false;
	  // here we loop over outliers as this will be faster than looping over the used hits
          for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
	    if ( rawHit.id() == outliers[ohit].first.id() ) { 
              is_outlier = true;                                    
              break; // break out of loop over outliers
            }
          }
          if (is_outlier == false) {
            used_hits.push_back(hit_list[ihit]);
            track->addToTrackerHits(used_hits.back());
            break; // break out of loop over rawObjects
          }          
        }
      } else {
        bool is_outlier = false;
        // here we loop over outliers as this will be faster than looping over the used hits
        for ( unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
          if ( trkHit == outliers[ohit].first ) { 
            is_outlier = true;                                    
            break; // break out of loop over outliers
          }
        }
        if (is_outlier == false) {
          used_hits.push_back(hit_list[ihit]);
          track->addToTrackerHits(used_hits.back());
        }          
      }
    }

    
//    ///////////////////////////////////////////////////////////////////////////
//    // We now need to find out at which point the fit is constrained 
//    // and therefore be able to provide well formed (pos. def.) cov. matrices
//    ///////////////////////////////////////////////////////////////////////////
//    
    
        
    ///////////////////////////////////////////////////////
    // first hit
    ///////////////////////////////////////////////////////
    
    edm4hep::TrackState* trkStateAtFirstHit = new edm4hep::TrackState() ;
    edm4hep::ConstTrackerHit firstHit = hits_in_fit.back().first;

    ///////////////////////////////////////////////////////
    // last hit
    ///////////////////////////////////////////////////////
    
    edm4hep::TrackState* trkStateAtLastHit = new edm4hep::TrackState() ;
    edm4hep::ConstTrackerHit lastHit = hits_in_fit.front().first;
          
    edm4hep::ConstTrackerHit last_constrained_hit = 0 ;     
    marlintrk->getTrackerHitAtPositiveNDF(last_constrained_hit);

    return_error = marlintrk->smooth(lastHit);
    
    if ( return_error != IMarlinTrack::success ) { 
      //streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: return_code for smoothing to " << lastHit << " = " << return_error << " NDF = " << ndf << std::endl;
      delete trkStateAtFirstHit;
      delete trkStateAtLastHit;
      return return_error ;
    }

    
    ///////////////////////////////////////////////////////
    // first create trackstate at IP
    ///////////////////////////////////////////////////////
    const edm4hep::Vector3d point; // nominal IP
    
    edm4hep::TrackState* trkStateIP = new edm4hep::TrackState() ;
      
    
    ///////////////////////////////////////////////////////
    // make sure that the track state can be propagated to the IP 
    ///////////////////////////////////////////////////////
    
    return_error = marlintrk->propagate(point, firstHit, *trkStateIP, chi2, ndf ) ;
    
    if ( return_error != IMarlinTrack::success ) { 
      //streamlog_out(DEBUG4) << "MarlinTrk::finaliseLCIOTrack: return_code for propagation = " << return_error << " NDF = " << ndf << std::endl;
      delete trkStateIP;
      delete trkStateAtFirstHit;
      delete trkStateAtLastHit;

      return return_error ;
    }
    
    trkStateIP->location = MarlinTrk::Location::AtIP;
    track->addToTrackStates(*trkStateIP);
    track->setChi2(chi2);
    track->setNdf(ndf);
    
              
    ///////////////////////////////////////////////////////
    // set the track states at the first and last hits 
    ///////////////////////////////////////////////////////    
    
    ///////////////////////////////////////////////////////
    // @ first hit
    ///////////////////////////////////////////////////////
    
    //streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> create TrackState AtFirstHit" << std::endl ;

    
    return_error = marlintrk->getTrackState(firstHit, *trkStateAtFirstHit, chi2, ndf ) ;
    
    if ( return_error == 0 ) {
      trkStateAtFirstHit->location = MarlinTrk::Location::AtFirstHit;
      track->addToTrackStates(*trkStateAtFirstHit);
    } else {
      //streamlog_out( WARNING ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at First Hit " << firstHit << std::endl ;
      delete trkStateAtFirstHit;
    }
    
    double r_first = firstHit.getPosition()[0]*firstHit.getPosition()[0] + firstHit.getPosition()[1]*firstHit.getPosition()[1];
    
    track->setRadiusOfInnermostHit(sqrt(r_first));
    
    if ( atLastHit == 0 && atCaloFace == 0 ) {
    
      ///////////////////////////////////////////////////////
      // @ last hit
      ///////////////////////////////////////////////////////  
      
      //streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> create TrackState AtLastHit : using trkhit " << last_constrained_hit << std::endl ;
      
      edm4hep::Vector3d last_hit_pos(lastHit.getPosition());
      
      return_error = marlintrk->propagate(last_hit_pos, last_constrained_hit, *trkStateAtLastHit, chi2, ndf);
            
//      return_error = marlintrk->getTrackState(lastHit, *trkStateAtLastHit, chi2, ndf ) ;
      
      if ( return_error == 0 ) {
        trkStateAtLastHit->location = MarlinTrk::Location::AtLastHit;
        track->addToTrackStates(*trkStateAtLastHit);
      } else {
        //streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Last Hit " << last_constrained_hit << std::endl ;
        delete trkStateAtLastHit;
      }
      
//      const EVENT::FloatVec& ma = trkStateAtLastHit->getCovMatrix();
//      
//      Matrix_Is_Positive_Definite( ma );

      ///////////////////////////////////////////////////////
      // set the track state at Calo Face 
      ///////////////////////////////////////////////////////
      
      edm4hep::TrackState trkStateCalo;
      bool tanL_is_positive = trkStateIP->tanLambda >0 ;
      
      return_error = createTrackStateAtCaloFace(marlintrk, &trkStateCalo, last_constrained_hit, tanL_is_positive);
//      return_error = createTrackStateAtCaloFace(marlintrk, trkStateCalo, lastHit, tanL_is_positive);
      
      if ( return_error == 0 ) {
        trkStateCalo.location = MarlinTrk::Location::AtCalorimeter;
        track->addToTrackStates(trkStateCalo);
      } else {
        //streamlog_out( WARNING ) << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Calo Face "  << std::endl ;
        //delete trkStateCalo;
      }
    } else {
      track->addToTrackStates(*atLastHit);
      track->addToTrackStates(*atCaloFace);
    }
    
    ///////////////////////////////////////////////////////
    // done
    ///////////////////////////////////////////////////////
    return return_error;
  }
  
  
  int createTrackStateAtCaloFace( IMarlinTrack* marlintrk, edm4hep::TrackState* trkStateCalo, edm4hep::ConstTrackerHit trkhit, bool tanL_is_positive ){
    
    //streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> createTrackStateAtCaloFace : using trkhit " << trkhit << " tanL_is_positive = " << tanL_is_positive << std::endl ;
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( marlintrk == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: IMarlinTrack == NULL ")  ) ;
    }
    
    if( trkStateCalo == 0 ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: TrackImpl == NULL ")  ) ;
    }
    
    if( !trkhit.isAvailable() ){
      throw EVENT::Exception( std::string("MarlinTrk::createTrackStateAtCaloFace: TrackHit == NULL ")  ) ;
    }
        
    int return_error = 0;
    
    double chi2 = -DBL_MAX;
    int ndf = 0;
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.reset() ;  // reset to 0
    
    encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::ECAL ;
    encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::barrel;
    encoder[lcio::ILDCellID0::layer]  = 0 ;
    
    int detElementID = 0;
    
    return_error = marlintrk->propagateToLayer(encoder.lowWord(), trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    
    if (return_error == IMarlinTrack::no_intersection ) { // try forward or backward
      if (tanL_is_positive) {
        encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::fwd;
      }
      else{
        encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::bwd;
      }
      return_error = marlintrk->propagateToLayer(encoder.lowWord(), trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    }
    
    if (return_error !=IMarlinTrack::success ) {
      //streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> createTrackStateAtCaloFace :  could not get TrackState at Calo Face: return_error = " << return_error << std::endl ;
    }
    
    
    return return_error;
    
    
  }
  
  void addHitNumbersToTrack(edm4hep::Track* track, std::vector<edm4hep::ConstTrackerHit>& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( track == 0 ){
      throw std::runtime_error( std::string("MarlinTrk::addHitsToTrack: TrackImpl == NULL ")  ) ;
    }
    std::map<int, int> hitNumbers; 
    
    hitNumbers[UTIL::ILDDetID::VXD] = 0;
    hitNumbers[UTIL::ILDDetID::SIT] = 0;
    hitNumbers[UTIL::ILDDetID::FTD] = 0;
    hitNumbers[UTIL::ILDDetID::TPC] = 0;
    hitNumbers[UTIL::ILDDetID::SET] = 0;
    hitNumbers[UTIL::ILDDetID::ETD] = 0;
    
    for(unsigned int j=0; j<hit_list.size(); ++j) {
      cellID_encoder.setValue(hit_list.at(j).getCellID()) ;
      int detID = cellID_encoder[UTIL::ILDCellID0::subdet];
      ++hitNumbers[detID];
      //std::cout << "debug: " << "Hit from Detector " << detID << std::endl;     
    }
    
    int offset = 2 ;
    if ( hits_in_fit == false ) { // all hit atributed by patrec
      offset = 1 ;
    }
    track->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::VXD]);
    track->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::FTD]);
    track->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::SIT]);
    track->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::TPC]);
    track->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::SET]);
    track->addToSubDetectorHitNumbers(hitNumbers[UTIL::ILDDetID::ETD]);
    //track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - offset ] = hitNumbers[lcio::ILDDetID::VXD];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - offset ] = hitNumbers[lcio::ILDDetID::FTD];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - offset ] = hitNumbers[lcio::ILDDetID::SIT];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - offset ] = hitNumbers[lcio::ILDDetID::TPC];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - offset ] = hitNumbers[lcio::ILDDetID::SET];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - offset ] = hitNumbers[lcio::ILDDetID::ETD];
  }
  
  void addHitNumbersToTrack(edm4hep::Track* track, std::vector<std::pair<edm4hep::ConstTrackerHit , double> >& hit_list, bool hits_in_fit, UTIL::BitField64& cellID_encoder){
    
    ///////////////////////////////////////////////////////
    // check inputs 
    ///////////////////////////////////////////////////////
    if( track == 0 ){
      throw std::runtime_error( std::string("MarlinTrk::addHitsToTrack: TrackImpl == NULL ")  ) ;
    }
    
    std::map<int, int> hitNumbers; 
    
    hitNumbers[lcio::ILDDetID::VXD] = 0;
    hitNumbers[lcio::ILDDetID::SIT] = 0;
    hitNumbers[lcio::ILDDetID::FTD] = 0;
    hitNumbers[lcio::ILDDetID::TPC] = 0;
    hitNumbers[lcio::ILDDetID::SET] = 0;
    hitNumbers[lcio::ILDDetID::ETD] = 0;
    
    for(unsigned int j=0; j<hit_list.size(); ++j) {
      cellID_encoder.setValue(hit_list.at(j).first.getCellID()) ;
      int detID = cellID_encoder[UTIL::ILDCellID0::subdet];
      ++hitNumbers[detID];
      //std::cout << "debug: Hit from Detector " << detID << std::endl;     
    }
    
    int offset = 2 ;
    if ( hits_in_fit == false ) { // all hit atributed by patrec
      offset = 1 ;
    }
    track->addToSubDetectorHitNumbers(hitNumbers[lcio::ILDDetID::VXD]);
    track->addToSubDetectorHitNumbers(hitNumbers[lcio::ILDDetID::FTD]);
    track->addToSubDetectorHitNumbers(hitNumbers[lcio::ILDDetID::SIT]);
    track->addToSubDetectorHitNumbers(hitNumbers[lcio::ILDDetID::TPC]);
    track->addToSubDetectorHitNumbers(hitNumbers[lcio::ILDDetID::SET]);
    track->addToSubDetectorHitNumbers(hitNumbers[lcio::ILDDetID::ETD]);
    //track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - offset ] = hitNumbers[lcio::ILDDetID::VXD];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - offset ] = hitNumbers[lcio::ILDDetID::FTD];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - offset ] = hitNumbers[lcio::ILDDetID::SIT];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - offset ] = hitNumbers[lcio::ILDDetID::TPC];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - offset ] = hitNumbers[lcio::ILDDetID::SET];
    //track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - offset ] = hitNumbers[lcio::ILDDetID::ETD];
    
  }
  
}
