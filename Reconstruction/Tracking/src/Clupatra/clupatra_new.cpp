#include "clupatra_new.h"
#include <set>
#include <vector>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

///---- GEAR ----
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/TPCModule.h"

#include "gear/BField.h"

#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackStateImpl.h"

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GearSvc/IGearSvc.h"
using namespace MarlinTrk ;

namespace lcio{
	// bits 0-15 are reserved !?
	const int  ILDTrackTypeBit::SEGMENT   = 16  ;
	const int  ILDTrackTypeBit::COMPOSITE = 17  ;
}

extern gear::GearMgr* gearMgr; // = _gear->getGearMgr();
extern RuntimeMap<clupatra_new::CluTrack*, MarlinTrk::IMarlinTrack*> MarTrkof;
extern RuntimeMap<edm4hep::Track, MarlinTrk::IMarlinTrack*> MarTrk_of_edm4hepTrack;
extern RuntimeMap<edm4hep::Track, clupatra_new::TrackInfoStruct*> TrackInfo_of_edm4hepTrack;

namespace clupatra_new{

	bool TrackZSort::operator()( edm4hep::Track l, edm4hep::Track r){
		return (  std::abs( TrackInfo_of_edm4hepTrack(l)->zAvg )   <   std::abs( TrackInfo_of_edm4hepTrack(r)->zAvg )  ) ;
	}

	void ComputeTrackerInfo::operator()( edm4hep::Track o )
	{
		edm4hep::Track lTrk  = o ;

		// compute z-extend of this track segment
		// const edm4hep::TrackerHitVec& hv = lTrk->getTrackerHits() ;

		float zMin =  1e99 ;
		float zMax = -1e99 ;
		float zAvg =  0. ;

		if( lTrk.trackerHits_size() >  1 ) {
			zMin = lTrk.getTrackerHits(            0  ).getPosition()[2] ;
			zMax = lTrk.getTrackerHits( lTrk.trackerHits_size() -1  ).getPosition()[2] ;
			zAvg = ( zMax + zMin ) / 2. ;
		}

		if( zMin > zMax ){ // swap
			float d = zMax ;
			zMax = zMin ;
			zMin = d  ;
		}

		if( ! TrackInfo_of_edm4hepTrack(lTrk) )
			TrackInfo_of_edm4hepTrack(lTrk) =  new TrackInfoStruct ;

		TrackInfo_of_edm4hepTrack(lTrk)->zMin = zMin ;
		TrackInfo_of_edm4hepTrack(lTrk)->zMax = zMax ;
		TrackInfo_of_edm4hepTrack(lTrk)->zAvg = zAvg ;

	}
	bool TrackCircleDistance:: operator()( nnclu::Element<ClupaPlcioTrack>* h0, nnclu::Element<ClupaPlcioTrack>* h1){

		edm4hep::Track trk0 = h0->first->edm4hepTrack ;
		edm4hep::Track trk1 = h1->first->edm4hepTrack ;

		const TrackInfoStruct* ti0 =  TrackInfo_of_edm4hepTrack(trk0);
		const TrackInfoStruct* ti1 =  TrackInfo_of_edm4hepTrack(trk1);


		/* FIXME debug Mingrui
		   streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id()
		   << "  (  ti0->zAvg > ti1->zAvg ) = " << (  ti0->zAvg > ti1->zAvg )
		   << std::endl ;
		   */


		// don't allow  overlaps in z !!!!
		float epsilon = 0. ;
		if(  ti0->zAvg > ti1->zAvg ){

			if( ti1->zMax > ( ti0->zMin + epsilon )  )
				return false ;

		} else {

			if( ti0->zMax > ( ti1->zMin + epsilon ) )
				return false ;

		}

		double tl0 = std::abs( getTanLambda(trk0) ) ;
		double tl1 = std::abs( getTanLambda(trk1) ) ;

		if( getTanLambda(trk0) * getTanLambda(trk1)   < 0. )
			return false ; // require the same sign


		double dtl = 2. * std::abs( tl0 - tl1 ) / ( tl0 + tl1 ) ;

		/* FIXME debug Mingrui
		   streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id()
		   << " dtl : " << dtl << "   std::abs( tl0 + tl1 ) = " <<  std::abs( tl0 + tl1 )
		   << " (  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  ) = " << (  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  )
		   << std::endl ;
		   */

		if(  dtl > 2.  * _dCut  &&  std::abs( tl0 + tl1 ) > 1.e-2  )
			return false ;
		// for very steep tracks (tanL < 0.001 ) tanL might differ largely for curlers due to multiple scattering

		double r0 = 1. / getOmega(trk0)   ;
		double r1 = 1. / getOmega(trk1)  ;

		double r0abs = std::abs( r0 ) ;
		double r1abs = std::abs( r1 ) ;


		double d0 = getD0(trk0) ;
		double d1 = getD0(trk1) ;

		// double z0 = trk0->getZ0() ;
		// double z1 = trk1->getZ0() ;

		double rIP = 20. ;

		// streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id()
		//  			      << " (  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  && std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     ) "
		//  			      << (  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  && std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     )
		//  			      << std::endl ;
		// // // don't merge tracks that come from an area of 20 mm around the IP
		// if(  std::abs( d0 ) < rIP &&  std::abs( z0 ) < rIP  &&
		//  	   std::abs( d1 ) < rIP &&  std::abs( z1 ) < rIP     )
		// 	return false ;


		edm4hep::TrackState ts0 =  getTrackStateAt(trk0, lcio::TrackState::AtFirstHit ) ;
		edm4hep::TrackState ts1 =  getTrackStateAt(trk1, lcio::TrackState::AtFirstHit ) ;

		/*
		   double z0 = ( ts0 ? ts0->referencePoint()[2]  : 0. ) ;
		   double z1 = ( ts1 ? ts1->referencePoint()[2]  : 0. ) ;
		   */
		// FIXME Mingrui since it is no pointer:
		double z0 = ts0.referencePoint[2];
		double z1 = ts1.referencePoint[2];


		/* FIXME debug Mingrui
		   streamlog_out( DEBUG2 ) << "TrackCircleDistance::operator() : " <<  trk0->id() << " <-> "  << trk1->id()
		   << " (  std::abs( z0 ) < rIP  &&  std::abs( z1 ) < rIP     ) "
		   << ( std::abs( z0 ) < rIP  &&  std::abs( z1 ) < rIP     )
		   << std::endl ;
		   */

		if(  std::abs( z0 ) < rIP  && std::abs( z1 ) < rIP     )
			return false ;



		double p0 = getPhi(trk0) ;
		double p1 = getPhi(trk1) ;

		double x0 = ( r0 - d0 ) * sin( p0 ) ;
		double x1 = ( r1 - d1 ) * sin( p1 ) ;

		double y0 = ( d0 - r0 ) * cos( p0 ) ;
		double y1 = ( d1 - r1 ) * cos( p1 ) ;

		double dr = 2. * std::abs( r0abs - r1abs )  / (r0abs + r1abs )  ;

		double distMS = sqrt ( ( x0 - x1 ) * ( x0 - x1 ) + ( y0 - y1 ) * ( y0 - y1 )  ) ;

		/* FIXME debug Mingrui
		   streamlog_out( DEBUG2 ) << "TrackCircleDistance:: operator() : " <<  trk0->id() << " <-> "  << trk1->id()
		   << "( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  ) "
		   << ( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  )
		   <<  " dr : " << dr
		   <<  " _dCut * std::abs( r0 ) " << _dCut * std::abs( r0 )
		   << " distMS :" << distMS
		   << std::endl ;
		   */

		//      return ( dr < DRMAX && distMS < _dCut * std::abs( r0 )  ) ;
		return ( dr < _dCut * std::abs( r0 )  &&  distMS < _dCut * std::abs( r0 )  ) ;

	}

	/** helper class to compute the chisquared of two points in rho and z coordinate */
	struct Chi2_RPhi_Z_Hit{
		//    double operator()( const TrackerHit* h, const gear::Vector3D& v1) {
		double operator()( const ClupaHit* h, const edm4hep::Vector3d& v1_t) {

			gear::Vector3D v1(v1_t[0], v1_t[1], v1_t[2]);
			//      gear::Vector3D v0( h->getPosition()[0] ,  h->getPosition()[1] ,  h->getPosition()[2] ) ;
			const gear::Vector3D& v0 = h->pos ;

			double sigsr =  h->edm4hepHit.getCovMatrix()[0] + h->edm4hepHit.getCovMatrix()[2]  ;
			double sigsz =  h->edm4hepHit.getCovMatrix()[5] ;
			// double sigsr =  0.01 ;
			// double sigsz =  0.1 ;

			double dPhi = std::abs(  v0.phi() - v1.phi() )  ;
			if( dPhi > M_PI )
				dPhi = 2.* M_PI - dPhi ;

			double dRPhi =  dPhi *  v0.rho() ;

			double dZ = v0.z() - v1.z() ;

			return  dRPhi * dRPhi / sigsr + dZ * dZ / sigsz  ;
		}
	};

	//-------------------------------------------------------------------------------

	struct CDot{
		CDot( int i , int j , double d) : I(i) , J(j) , Dot( d ) {}
		int I ;
		int J ;
		double Dot ;
	};

	bool sort_CDot( const CDot& c0 ,const CDot& c1 ){
		return c0.Dot > c1.Dot ;
	}

	//-------------------------------------------------------------------------------


	int addHitsAndFilter( CluTrack* clu, HitListVector& hLV , double dChi2Max, double chi2Cut, unsigned maxStep, ZIndex& zIndex, bool backward,
			MarlinTrk::IMarlinTrkSystem* trkSys ) {


		int nHitsAdded = 0 ;

		const double bfield = gearMgr->getBField().at( gear::Vector3D(0.,0.0,0.) ).z() ;

		// Support for more than one module
		static const gear::TPCParameters*  gearTPC = &(gearMgr->getTPCParameters());
		// The ternary operator is used to make the trick with the static variable which
		// is supposed to be calculated only once, also for performance reason
		static const int maxTPCLayerID  = (gearMgr->getDetectorName() == "LPTPC" ) ?
			gearTPC->getModule(0).getNRows() + gearTPC->getModule(2).getNRows() + gearTPC->getModule(5).getNRows() - 1 : // LCTPC
			gearTPC->getModule(0).getNRows() - 1 ; // ILD

		clu->sort( LayerSortIn() ) ;

		int layer =  ( backward ?  clu->front()->first->layer : clu->back()->first->layer   ) ;


		// std::cout <<  " ======================  addHitsAndFilter():  - layer " << layer << "  backward: " << backward << std::endl;
		// std::cerr <<  " ======================  addHitsAndFilter():  - layer " << layer << "  backward: " << backward << std::endl;


		if( layer <= 0  || layer >=  maxTPCLayerID   )
			return  nHitsAdded ;


		Chi2_RPhi_Z_Hit ch2rzh ;

		IMarlinTrack* trk =  MarTrkof(clu);
// {
// edm4hep::TrackState ts;
// double chi2; int ndf;
// trk->getTrackState( ts, chi2,  ndf ) ;
// std::cout << "Omega of this part is: " << ts.omega << " " << chi2 << " " << ndf << std::endl;
// }

		if(  trk == 0 ){
			// streamlog_out( DEBUG3 ) <<  "  addHitsAndFilter called with null pointer to MarlinTrk  - won't do anything " << std::endl ;
			return  nHitsAdded;
		}

		unsigned step = 0 ;

		UTIL::BitField64 encoder( UTIL::ILDCellID0::encoder_string ) ;
		encoder[UTIL::ILDCellID0::subdet] = UTIL::ILDDetID::TPC ;

		edm4hep::ConstTrackerHit firstHit; // =  0 ;
                firstHit.unlink();

		IMarlinTrack* bwTrk = 0 ;

		if( trkSys && backward  ) { //==================== only active if called with _trkSystem pointer ============================
                        // std::cout << "It is backward" << std::endl;

			// need to go back in cluster until 4th hit from the start
			CluTrack::iterator it =  clu->begin() , end = clu->end() ;
			int i=0 ;
			for(     ; it++ != end && i< 3 ;  ++i  ) ;

			if( it == end ) --it ;

			// std::cout <<  " ---- addHitsAndFilter : will smooth back to " << i <<"th  hit - size of clu " << clu->size() << std::endl ;

			if( !  (*it)->first ){

				// std::cout << " ---- addHitsAndFilter : null pointer at i-th hit : i = " << i << std::endl ;
				return nHitsAdded ;
			}


			firstHit = (*it)->first->edm4hepHit ;

			int smoothed  = trk->smooth( firstHit ) ;

			double chi2 ;
			int ndf  ;
			edm4hep::TrackState ts ;
			trk->getTrackState( firstHit , ts, chi2,  ndf ) ;

			// std::cout <<  "  -- addHitsAndFilter(): smoothed track segment : " <<  MarlinTrk::errorCode( smoothed )
			//  <<  " using track state : " <<   ts
			//  <<  " ---- chi2: " << chi2
			//  <<  "  ndf " << ndf
			//  << std::endl ;



			//      std::auto_ptr<MarlinTrk::IMarlinTrack> mTrk( _trksystem->createTrack()  ) ;
			bwTrk = trkSys->createTrack()  ;

			//need to add a dummy hit to the track
			bwTrk->addHit(  firstHit  ) ; // use the hit we smoothed back to

			// note: backward here is forward in the MarlinTrk sense, ie. in direction of energy loss
			bwTrk->initialise( ts ,  bfield ,  MarlinTrk::IMarlinTrack::forward ) ;

		}


		IMarlinTrack* theTrk  = ( bwTrk ? bwTrk : trk )  ;


		while( step < maxStep + 1 ) {

			layer += ( backward ?  +1 : -1 )  ;

			if( layer < 0  || layer >  maxTPCLayerID   )
				break ;


			encoder[ UTIL::ILDCellID0::layer ]  = layer ;

			// gear::Vector3D xv ;
			edm4hep::Vector3d xv;

			bool hitAdded = false ;

			int layerID = encoder.lowWord() ;
			int elementID = 0 ;

			//      int mode = ( backward ? IMarlinTrack::modeBackward : IMarlinTrack::modeForward  )  ;
			//      int mode = ( backward ? IMarlinTrack::modeForward : IMarlinTrack::modeClosest )  ;
			int mode =  IMarlinTrack::modeClosest ;


			int intersects = -1  ;

			if( firstHit.isAvailable() )  { // FIXME Mingrui:: Important!! need to check whether firstHit is OK or not.
// std::cout << "available" << std::endl;
// std::cout << layerID << std::endl;
// std::cout << firstHit << std::endl;
// std::cout << xv << std::endl;
// std::cout << elementID << std::endl;
// std::cout << mode << std::endl;
			intersects  = theTrk->intersectionWithLayer( layerID, firstHit, xv, elementID , mode )   ;

			} else {

			 	intersects  = theTrk->intersectionWithLayer( layerID, xv, elementID , mode )  ;
			}



			   // std::cout <<  "  -- addHitsAndFilter(): looked for intersection - "
			   // <<  "  Step : " << step
			   // <<  "  at layer: "   << layer
			   // <<  "   intersects: " << MarlinTrk::errorCode( intersects )
			   // <<  "  next xing point : " <<  xv  ;


			if( intersects == IMarlinTrack::success ) { // found a crossing point

				int zIndCP = zIndex.index( xv[2] ) ;

				HitList& hLL = hLV.at( layer ) ;

				double ch2Min = 1.e99 ;
				Hit* bestHit = 0 ;

				// streamlog_out( DEBUG3 ) <<  "      -- number of hits on layer " << layer << " : " << hLL.size() << std::endl ;

				for( HitList::const_iterator ih = hLL.begin(), end = hLL.end() ; ih != end ; ++ih ){

					// if the z indices differ by more than one we can continue
					if( nnclu::notInRange<-1,1>(  (*ih)->first->zIndex - zIndCP ) )
						continue ;

					double ch2 = ch2rzh( (*ih)->first , xv )  ;

					if( ch2 < ch2Min ){

						ch2Min = ch2 ;
						bestHit = (*ih) ;
					}
				}//-------------------------------------------------------------------

				if( bestHit != 0 ){

					// streamlog_out( DEBUG3 ) <<   " ************ bestHit "  << bestHit
					// <<   " pos : " <<   (bestHit ? bestHit->first->pos :  gear::Vector3D() )
					// <<   " chi2: " <<  ch2Min
					// <<   " chi2Cut: " <<  chi2Cut <<   std::endl ;

					const gear::Vector3D&  hPos = bestHit->first->pos  ;

					if( ch2Min  < chi2Cut ) {

						double deltaChi = 0. ;

						edm4hep::ConstTrackerHit ht = bestHit->first->edm4hepHit;
						int addHit =  theTrk->addAndFit(ht , deltaChi, dChi2Max )  ;



						// streamlog_out( DEBUG3 ) <<   " *****       assigning left over hit : " << errorCode( addHit )  //<< hPos << " <-> " << xv
						/*
						   <<   " dist: " <<  (  hPos - xv ).r()
						   <<   " chi2: " <<  ch2Min
						   <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->edm4hepHit->getCovMatrix()[0]
						   + bestHit->first->edm4hepHit->getCovMatrix()[2] )
						   <<	 "  z= " <<  sqrt( bestHit->first->edm4hepHit->getCovMatrix()[5] )
						   <<   "  ----- deltaChi = " << deltaChi
						   << std::endl ;
						   */


						if( addHit  == IMarlinTrack::success ) {

							hitAdded = true ;

							hLL.remove(  bestHit ) ;
							clu->addElement( bestHit ) ;

							// firstHit = 0 ; // after we added a hit, the next intersection search should use this last hit...
							firstHit.unlink();

							++nHitsAdded ;

							// streamlog_out( DEBUG ) <<   " ---- track state filtered with new hit ! ------- " << std::endl ;
						}
					} // chi2Cut
				} // bestHit


			} else {  // no intersection found

				// we stop searching here
				step = maxStep + 1 ;
			}

			// reset the step or increment it
			step = (  hitAdded  ?  1 :  step + 1 ) ;

		} // while step < maxStep


		if( bwTrk ) delete bwTrk ;

		return nHitsAdded ;

	}

	//------------------------------------------------------------------------------------------------------------

	bool addHitAndFilter( int detectorID, int layer, CluTrack* clu, HitListVector& hLV , double dChi2Max, double chi2Cut) {

		Chi2_RPhi_Z_Hit ch2rzh ;

		IMarlinTrack* trk =  MarTrkof(clu) ;

		UTIL::BitField64 encoder( UTIL::ILDCellID0::encoder_string ) ;

		encoder[ UTIL::ILDCellID0::subdet ] = detectorID ;
		encoder[ UTIL::ILDCellID0::layer  ] = layer ;

		int layerID = encoder.lowWord() ;

		// gear::Vector3D xv ;
		edm4hep::Vector3d xv ;

		int elementID = -1 ;

		int intersects = trk->intersectionWithLayer( layerID, xv, elementID , IMarlinTrack::modeClosest  ) ;

		//  IMarlinTrack::modeBackward , IMarlinTrack::modeForward

		/*
		   streamlog_out( DEBUG2 ) <<  "  ============ addHitAndFilter(): looked for intersection - "
		   <<  "  detector : " <<  detectorID
		   <<  "  at layer: "   << layer
		   <<  "  intersects: " << MarlinTrk::errorCode( intersects )
		   <<  "  next xing point : " <<  xv  ;
		   */

		bool hitAdded = false ;

		if( intersects == IMarlinTrack::success ) { // found a crossing point

			//FIXME: this is copied from ClupatraProcessor
			static const int ZBins = 160 ;
			ZIndex zIndex( -2750. , 2750. ,ZBins  ) ;
			int zIndCP = zIndex.index( xv[2] ) ;

			HitList& hLL = hLV.at( layer ) ;

			double ch2Min = 1.e99 ;
			Hit* bestHit = 0 ;

			for( HitList::const_iterator ih = hLL.begin(), end = hLL.end() ; ih != end ; ++ih ){

				// if the z indices differ by more than one we can continue
				if( nnclu::notInRange<-1,1>(  (*ih)->first->zIndex - zIndCP ) )
					continue ;

				double ch2 = ch2rzh( (*ih)->first , xv )  ;

				if( ch2 < ch2Min ){

					ch2Min = ch2 ;
					bestHit = (*ih) ;
				}
			}//-------------------------------------------------------------------


			// streamlog_out( DEBUG2 ) <<   " ************ bestHit "  << bestHit
			// <<   " pos : " <<   (bestHit ? bestHit->first->pos :  gear::Vector3D() )
			// <<   " chi2: " <<  ch2Min
			// <<   " chi2Cut: " <<  chi2Cut <<   std::endl ;

			if( bestHit != 0 ){

				const gear::Vector3D&  hPos = bestHit->first->pos  ;

				if( ch2Min  < chi2Cut ) {

					double deltaChi = 0. ;
					edm4hep::ConstTrackerHit bh = bestHit->first->edm4hepHit;
					int addHit = trk->addAndFit( bh, deltaChi, dChi2Max ) ;



					// streamlog_out( DEBUG2 ) <<   " *****       assigning left over hit : " << errorCode( addHit )  //<< hPos << " <-> " << xv
					/*
					   <<   " dist: " <<  (  hPos - xv ).r()
					   <<   " chi2: " <<  ch2Min
					   <<   "  hit errors :  rphi=" <<  sqrt( bestHit->first->edm4hepHit->getCovMatrix()[0]
					   + bestHit->first->edm4hepHit->getCovMatrix()[2] )
					   <<	 "  z= " <<  sqrt( bestHit->first->edm4hepHit->getCovMatrix()[5] )
					   << std::endl ;
					   */





					if( addHit  == IMarlinTrack::success ) {

						hitAdded = true ;

						hLL.remove(  bestHit ) ;
						clu->addElement( bestHit ) ;


						// streamlog_out( DEBUG2 ) <<   " ---- track state filtered with new hit ! ------- " << std::endl ;
					}
				} // chi2Cut
			} // bestHit

		} // intersection found

		return hitAdded ;
	}


	//------------------------------------------------------------------------------------------------------------

	void getHitMultiplicities( CluTrack* clu, std::vector<int>& mult ){

		// int static maxTPCLayerID = marlin::Global::GEAR->getTPCParameters().getPadLayout().getNRows() - 1 ;
		// HitListVector hLV( maxTPCLayerID )  ;
		// addToHitListVector( clu->begin(), clu->end(),  hLV ) ;

		std::map<int, unsigned > hm ;

		for( CluTrack::iterator it=clu->begin(), end =clu->end() ;   it != end ; ++ it ){

			++hm[ (*it)->first->layer ] ;
		}

		unsigned maxN = mult.size() - 1 ;

		//    for( unsigned i=0 ; i < hLV.size() ; ++i ){
		//      unsigned m =  hLV[i].size() ;

		for( std::map<int, unsigned>::iterator it= hm.begin() , end = hm.end() ; it != end ; ++ it ){

			unsigned m = ( it->second < maxN ?  it->second  :  maxN   ) ;

			if( m == 0 )
				continue ;

			++mult[m] ;

			++mult[0] ;
		}
	}

	//-------------------------------------------------------------------------------------------------------------------------
	void split_multiplicity( Clusterer::cluster_list& cluList, int layerWithMultiplicity , int N) {

		for( Clusterer::cluster_list::iterator it= cluList.begin(), end= cluList.end() ; it != end ; ++it ){

			CluTrack* clu = *it ;

			std::vector<int> mult( N + 2 ) ;
			// get hit multiplicities up to N ( N+1 means N+1 or higher )
			getHitMultiplicities( clu , mult ) ;


			// streamlog_out(  DEBUG2 ) << " **** split_multiplicity -  hit multiplicities: \n" ;

			/*
			   for( unsigned i=0,n=mult.size() ; i<n ; ++i) {
			   streamlog_out(  DEBUG2 ) << "     m["<<i<<"] = " <<  mult[i] << "\n"  ;
			   }
			   */


			if(  mult[1] >= layerWithMultiplicity )
				continue ;

			bool split_cluster = false  ;

			for( int m=2 ; m<N ; ++m){

				if( m == 2 && mult[m] >= layerWithMultiplicity ){

					// streamlog_out(  DEBUG3 ) << " **** split_multiplicity - create_two_clusters \n" ;

					create_two_clusters( *clu , cluList  ) ;

					split_cluster = true  ;
				}
				else if( m == 3 && mult[m] >= layerWithMultiplicity ){

					// streamlog_out(  DEBUG3 ) << " **** split_multiplicity - create_three_clusters \n" ;

					create_three_clusters( *clu , cluList  ) ;

					split_cluster = true  ;
				}
				else if(  mult[m] >= layerWithMultiplicity ){

					// streamlog_out(  DEBUG3 ) << " **** split_multiplicity - create_n_clusters \n" ;

					create_n_clusters( *clu ,cluList , m  ) ;

					split_cluster = true  ;
				}

			}
			if( split_cluster ){

				clu->clear() ;

				it = cluList.erase( it ) ;
				--it ; // erase returns iterator to next element
			}


		}
	}

	//-------------------------------------------------------------------------------------------------------------------------

	void create_n_clusters( Clusterer::cluster_type& hV, Clusterer::cluster_list& cluVec , unsigned n ) {

		if( n < 4 ){

			// streamlog_out( ERROR ) <<  " create_n_clusters called for n = " << n << std::endl ;
			return ;
		}

		// copy of the algorithm for creating three clusters (minimize angle between hits as seen from IP)

		hV.freeElements() ;

		// Support for more than one module
		//auto _gear = service<IGearSvc>("GearSvc");

		static const gear::TPCParameters*  gearTPC = &(gearMgr->getTPCParameters());
		// The ternary operator is used to make the trick with the static variable which
		// is supposed to be calculated only once, also for performance reason
		static const int tpcNRow =
			(gearMgr->getDetectorName() == "LPTPC" ) ?
			gearTPC->getModule(0).getNRows() + gearTPC->getModule(2).getNRows() + gearTPC->getModule(5).getNRows() : // LCTPC
			gearTPC->getModule(0).getNRows() ; // ILD

		HitListVector hitsInLayer( tpcNRow )  ;
		addToHitListVector(  hV.begin(), hV.end(), hitsInLayer ) ;

		std::vector< CluTrack*> clu(n)  ;

		std::vector< gear::Vector3D > lastp(n) ;

		for(unsigned i=0; i<n; ++i){

			clu[i] = new CluTrack ;

			cluVec.push_back( clu[i] ) ;
		}

		for( int l=tpcNRow-1 ; l >= 0 ; --l){

			HitList& hL = hitsInLayer[ l ] ;

			// streamlog_out(  DEBUG ) << " create_n_clusters  --- layer " << l  <<  " size: " << hL.size() << std::endl ;

			if( hL.size() != n ){  // ignore layers with different hit numbers

				continue ;
			}

			HitList::iterator iH = hL.begin() ;

			std::vector<Hit*> h(n) ;
			std::vector< gear::Vector3D>  p(n) ;


			// streamlog_out(  DEBUG2 ) << " create_n_clusters  ---  layer " << l << std::endl ;

			for(unsigned i=0; i<n; ++i){

				h[i] = *iH++ ;

				p[i] = h[i]->first->pos ;

				// streamlog_out(  DEBUG2 ) << "              h"<< i << ": " << p[i] << std::endl ;
			}

			// streamlog_out(  DEBUG2 ) << " ---------------------------- " << std::endl ;


			if( clu[0]->empty() ){ // first hit tuple

				// streamlog_out(  DEBUG ) << " create_n_clusters  --- initialize clusters " << std::endl ;

				for(unsigned i=0; i<n; ++i){

					clu[i]->addElement( h[i] ) ;

					lastp[i] = ( 1. / p[i].r() ) * p[i] ;
				}

				continue ;   //----------  that's all for the first layer w/ hits
			}


			// unit vector in direction of current hits
			std::vector< gear::Vector3D > pu(n) ;

			for(unsigned i=0; i<n; ++i){

				pu[i] = ( 1. / p[i].r() ) * p[i] ;
			}

			std::list< CDot > cDot ;

			// create list of dot products between last hit in cluster and current hit
			// cos( angle between  hits as seen from IP )
			for( unsigned i = 0 ; i < n ; ++i ){
				for( unsigned j = 0 ; j < n ; ++j ){

					cDot.push_back( CDot( i , j , lastp[i].dot( pu[ j ] ) )) ;
				}
			}

			// sort dot products in descending order
			cDot.sort( sort_CDot ) ;

			// assign hits to clusters with largest dot product ( smallest angle )

			std::set<int> cluIdx ;
			std::set<int> hitIdx ;

			for( std::list<CDot>::iterator it = cDot.begin() ; it != cDot.end() ; ++ it ) {
				/*
				   streamlog_out(  DEBUG2 ) << " I : " << it->I
				   << " J : " << it->J
				   << " d : " << it->Dot
				   << std::endl ;
				   */
				unsigned i =  it->I ;
				unsigned j =  it->J ;

				// ignore clusters or hits that hvae already been assigned
				if( cluIdx.find( i ) == cluIdx.end()  && hitIdx.find( j ) == hitIdx.end() ){

					cluIdx.insert( i ) ;
					hitIdx.insert( j ) ;

					clu[ i ]->addElement( h[ j ] ) ;

					// cluDir[i] = - ( p[j] - lastp[i] ) ;
					// cluDir[i] =  (  1. / cluDir[i].r() ) * cluDir[i]  ;

					lastp[i] =  p[j] ;

					/*
					   streamlog_out(  DEBUG2 ) << " **** adding to cluster : " << it->I
					   << " hit  : " << it->J
					   << " d : " << it->Dot
					   << std::endl ;
					   */
				}

			}
		}


		// streamlog_out(  DEBUG ) << " create_three_clusters  --- clu[0] " << clu[0]->size()
		// <<  " clu[1] " << clu[1]->size()
		// <<  " clu[2] " << clu[2]->size()
		// << std::endl ;




		return ;
	}


	//======================================================================================================================

	void create_three_clusters( Clusterer::cluster_type& hV, Clusterer::cluster_list& cluVec ) {

		hV.freeElements() ;

		// auto _gear = service<IGearSvc>("GearSvc");
		// Support for more than one module
		static const gear::TPCParameters*  gearTPC = &(gearMgr->getTPCParameters());
		// The ternary operator is used to make the trick with the static variable which
		// is supposed to be calculated only once, also for performance reason
		static const int tpcNRow =
			(gearMgr->getDetectorName() == "LPTPC" ) ?
			gearTPC->getModule(0).getNRows() + gearTPC->getModule(2).getNRows() + gearTPC->getModule(5).getNRows() : // LCTPC
			gearTPC->getModule(0).getNRows() ; // ILD

		HitListVector hitsInLayer( tpcNRow )  ;
		addToHitListVector(  hV.begin(), hV.end(), hitsInLayer ) ;

		CluTrack* clu[3] ;

		gear::Vector3D lastp[3] ;
		//    gear::Vector3D cluDir[3] ;

		clu[0] = new CluTrack ;
		clu[1] = new CluTrack ;
		clu[2] = new CluTrack ;

		cluVec.push_back( clu[0] ) ;
		cluVec.push_back( clu[1] ) ;
		cluVec.push_back( clu[2] ) ;


		for( int l=tpcNRow-1 ; l >= 0 ; --l){

			HitList& hL = hitsInLayer[ l ] ;

			// streamlog_out(  DEBUG ) << " create_three_clusters  --- layer " << l  <<  " size: " << hL.size() << std::endl ;

			if( hL.size() != 3 ){

				continue ;
			}

			HitList::iterator iH = hL.begin() ;

			Hit* h[3] ;
			h[0] = *iH++ ;
			h[1] = *iH++   ;
			h[2] = *iH   ;

			gear::Vector3D p[3] ;
			p[0] = h[0]->first->pos ;
			p[1] = h[1]->first->pos ;
			p[2] = h[2]->first->pos ;


			/*
			   streamlog_out(  DEBUG ) << " create_three_clusters  ---  layer " << l
			   << "  h0 : " << p[0]
			   << "  h1 : " << p[1]
			   << "  h2 : " << p[2] ;
			   */
			//			       << std::endl ;

			if( clu[0]->empty() ){ // first hit triplet

				// streamlog_out(  DEBUG ) << " create_three_clusters  --- initialize clusters " << std::endl ;

				clu[0]->addElement( h[0] ) ;
				clu[1]->addElement( h[1] ) ;
				clu[2]->addElement( h[2] ) ;

				// lastp[0] =  p[0] ;
				// lastp[1] =  p[1] ;
				// lastp[2] =  p[2] ;

				lastp[0] = ( 1. / p[0].r() ) * p[0] ;
				lastp[1] = ( 1. / p[1].r() ) * p[1] ;
				lastp[2] = ( 1. / p[2].r() ) * p[2] ;

				continue ;
			}


			// unit vector in direction of current hits
			gear::Vector3D pu[3] ;
			pu[0] = ( 1. / p[0].r() ) * p[0] ;
			pu[1] = ( 1. / p[1].r() ) * p[1] ;
			pu[2] = ( 1. / p[2].r() ) * p[2] ;


			std::list< CDot > cDot ;

			// create list of dot products between last hit in cluster and current hit
			// cos( angle between  hits as seen from IP )
			for( int i = 0 ; i < 3 ; ++i ){
				for( int j = 0 ; j < 3 ; ++j ){

					// // unit vector in direction of last hit and
					// gear::Vector3D pu = - ( p[0] - lastp[0] ) ;
					// pu[1] = - ( p[1] - lastp[1] ) ;
					// pu[2] = - ( p[2] - lastp[2] ) ;
					// pu[0] =  (  1. / pu[0].r() ) * pu[0]  ;
					// pu[1] =  (  1. / pu[1].r() ) * pu[1]  ;
					// pu[2] =  (  1. / pu[2].r() ) * pu[2]  ;
					// cDot.push_back( CDot( i , j , cluDir[i].dot( pu[ j ] ) )) ;

					cDot.push_back( CDot( i , j , lastp[i].dot( pu[ j ] ) )) ;

				}
			}

			// sort dot products in descending order
			cDot.sort( sort_CDot ) ;

			// assign hits to clusters with largest dot product ( smallest angle )

			std::set<int> cluIdx ;
			std::set<int> hitIdx ;

			for( std::list<CDot>::iterator it = cDot.begin() ; it != cDot.end() ; ++ it ) {

				/*
				   streamlog_out(  DEBUG ) << " I : " << it->I
				   << " J : " << it->J
				   << " d : " << it->Dot
				   << std::endl ;
				   */

				int i =  it->I ;
				int j =  it->J ;

				// ignore clusters or hits that hvae already been assigned
				if( cluIdx.find( i ) == cluIdx.end()  && hitIdx.find( j ) == hitIdx.end() ){

					cluIdx.insert( i ) ;
					hitIdx.insert( j ) ;

					clu[ i ]->addElement( h[ j ] ) ;

					// cluDir[i] = - ( p[j] - lastp[i] ) ;
					// cluDir[i] =  (  1. / cluDir[i].r() ) * cluDir[i]  ;

					lastp[i] =  p[j] ;

					/*
					   streamlog_out(  DEBUG ) << " **** adding to cluster : " << it->I
					   << " hit  : " << it->J
					   << " d : " << it->Dot
					   << std::endl ;
					   */
				}

			}

			// lastp[0] =  pu[0] ;
			// lastp[1] =  pu[1] ;
			// lastp[2] =  pu[2] ;
			// lastp[0] =  p[0] ;
			// lastp[1] =  p[1] ;
			// lastp[2] =  p[2] ;


		}


		// streamlog_out(  DEBUG ) << " create_three_clusters  --- clu[0] " << clu[0]->size()
		//<<  " clu[1] " << clu[1]->size()
		//<<  " clu[2] " << clu[2]->size()
		//<< std::endl ;




		return ;
	}
	//-----------------------------------------------------------------

	void create_two_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ) {


		clu.freeElements() ;

		// streamlog_out(  DEBUG ) << " create_two_clusters  --- called ! - size :  " << clu.size()  << std::endl ;
		// auto _gear = service<IGearSvc>("GearSvc");


		// Support for more than one module
		static const gear::TPCParameters*  gearTPC = &(gearMgr->getTPCParameters());
		// The ternary operator is used to make the trick with the static variable which
		// is supposed to be calculated only once, also for performance reason
		static const int tpcNRow =
			(gearMgr->getDetectorName() == "LPTPC" ) ?
			gearTPC->getModule(0).getNRows() + gearTPC->getModule(2).getNRows() + gearTPC->getModule(5).getNRows() : // LCTPC
			gearTPC->getModule(0).getNRows() ; // ILD

		HitListVector hitsInLayer( tpcNRow )  ;

		addToHitListVector(  clu.begin(), clu.end(), hitsInLayer ) ;


		CluTrack* clu0 = new CluTrack ;
		CluTrack* clu1 = new CluTrack ;

		cluVec.push_back( clu0 ) ;
		cluVec.push_back( clu1 ) ;


		gear::Vector3D lastDiffVec(0.,0.,0.) ;

		for( int l=tpcNRow-1 ; l >= 0 ; --l){

			HitList& hL = hitsInLayer[ l ] ;

			// streamlog_out(  DEBUG ) << " create_two_clusters  --- layer " << l  <<  " size: " << hL.size() << std::endl ;

			if( hL.size() != 2 ){

				// if there is not exactly two hits in the layer, we simply copy the unassigned hits to the leftover hit vector
				// std::copy( hL.begin(),  hL.end() , std::back_inserter( hV ) ) ;
				continue ;
			}

			HitList::iterator iH = hL.begin() ;

			Hit* h0 = *iH++ ;
			Hit* h1 = *iH   ;

			gear::Vector3D& p0 = h0->first->pos ;
			gear::Vector3D& p1 = h1->first->pos ;


			// streamlog_out(  DEBUG ) << " create_two_clusters  ---  layer " << l
			//<< "  h0 : " << p0
			// << "  h1 : " << p1 ;
			//			       << std::endl ;

			if( clu0->empty() ){ // first hit pair

				// streamlog_out(  DEBUG ) << " create_two_clusters  --- initialize clusters " << std::endl ;

				clu0->addElement( h0 ) ;
				clu1->addElement( h1 ) ;

				lastDiffVec = p1 - p0 ;

				continue ;
			}

			gear::Vector3D d = p1 - p0 ;

			float s0 =  ( lastDiffVec + d ).r() ;
			float s1 =  ( lastDiffVec - d ).r() ;

			if( s0 > s1 ){  // same orientation, i.e. h0 in this layer belongs to h0 in first layer

				// streamlog_out(  DEBUG ) << " create_two_clusters  ---   same orientation " << std::endl ;
				clu0->addElement( h0 ) ;
				clu1->addElement( h1 ) ;

			} else{                // oposite orientation, i.e. h1 in this layer belongs to h0 in first layer

				// streamlog_out(  DEBUG2 ) << " create_two_clusters  ---  oposite orientation " << std::endl ;
				clu0->addElement( h1 ) ;
				clu1->addElement( h0 ) ;
			}

		}

		// clu.freeHits() ;

		// streamlog_out(  DEBUG1 ) << " create_two_clusters  --- clu0 " << clu0->size()
		//  <<  " clu1 " << clu1->size() << std::endl ;

		// return std::make_pair( clu0 , clu1 ) ;

		return ;
	}

	//-------------------------------------------------------------------------------------------------------------------------


	// /** Find the nearest hits in the previous and next layers - (at most maxStep layers appart - NOT YET).
	//  */
	// void  findNearestHits( CluTrack& clu,  int maxLayerID, int maxStep){

	//   HitListVector hitsInLayer( maxLayerID )  ;
	//   addToHitListVector(  clu.begin(), clu.end(), hitsInLayer ) ;

	//   for(CluTrack::iterator it= clu.begin() ; it != clu.end() ; ++it ){

	//     Hit* hit0 = *it ;
	//     gear::Vector3D& pos0 =  hit0->first->pos ;

	//     int layer = hit0->first->layerID ;

	//     int l=0 ;
	//     if( (l=layer+1) <  maxLayerID ) {

	// 	HitList& hL = hitsInLayer[ l ] ;

	// 	double minDist2 = 1.e99 ;
	// 	Hit*   bestHit = 0 ;

	// 	for( HitList::iterator iH = hL.begin() ; iH != hL.end() ; ++iH ) {

	// 	  Hit* hit1 = *iH ;
	// 	  gear::Vector3D& pos1 = hit1->first->pos ;

	// 	  double dist2 = ( pos0 - pos1 ).r2() ;

	// 	  if( dist2 < minDist2 ){
	// 	    minDist2 = dist2 ;
	// 	    bestHit = hit1 ;
	// 	  }
	// 	}
	// 	if( bestHit ){

	// 	  hit0->first->nNNHit =  bestHit->first ;
	// 	  hit0->first->nDist  =  minDist2 ;
	// 	}
	//     }

	//     if( (l=layer-1) > 0  ) {

	// 	HitList& hL = hitsInLayer[ l ] ;

	// 	double minDist2 = 1.e99 ;
	// 	Hit* bestHit = 0 ;

	// 	for( HitList::iterator iH = hL.begin() ; iH != hL.end() ; ++iH ) {

	// 	  Hit* hit1 = *iH ;
	// 	  gear::Vector3D&  pos1 =  hit1->first->pos ;

	// 	  double dist2 = ( pos0 - pos1 ).r2() ;

	// 	  if( dist2 < minDist2 ){
	// 	    minDist2 = dist2 ;
	// 	    bestHit = hit1 ;
	// 	  }
	// 	}
	// 	if( bestHit ){

	// 	  hit0->first->pNNHit =  bestHit->first ;
	// 	  hit0->first->pDist  =  minDist2 ;
	// 	}
	//     }

	//   }
	// }


	//-------------------------------------------------------------------------------------------------------------------------

	MarlinTrk::IMarlinTrack* IMarlinTrkFitter::operator() (CluTrack* clu) {

		bool isFirstFit = true ;
		double maxChi2  =  _maxChi2Increment   ;

start:

		MarlinTrk::IMarlinTrack* trk = _ts->createTrack();

		//if( clu->empty()  ){
		if( clu->size() < 3  ){

			// streamlog_out( ERROR ) << " IMarlinTrkFitter::operator() : cannot fit cluster track with less than 3 hits ! " << std::endl ;

			return trk ;
		}

		MarTrkof(clu) = trk ;

		clu->sort( LayerSortOut() ) ;

		// need to reverse the order for incomming track segments (curlers)
		// assume particle comes from IP
		Hit* hf = clu->front() ;
		Hit* hb = clu->back() ;

		bool reverse_order =   ( std::abs( hf->first->pos.z() ) > std::abs( hb->first->pos.z()) + 3. ) ;

		unsigned nHit = 0 ;

		if( reverse_order ){
		  //std::cout << "It is true order" << std::endl;
		  for( CluTrack::reverse_iterator it=clu->rbegin() ; it != clu->rend() ; ++it){
		    edm4hep::ConstTrackerHit ph = (*it)->first->edm4hepHit;
		    trk->addHit(ph) ;
		    ++nHit ;
		    //std::cout  <<  "   hit  added  " <<  (*it)->first->edm4hepHit   << std::endl ;
		  }
		  
		  trk->initialise( MarlinTrk::IMarlinTrack::forward ) ;

		} else {
		  //std::cout << "It is reverse order" << std::endl;
		  for( CluTrack::iterator it=clu->begin() ; it != clu->end() ; ++it){
		    edm4hep::ConstTrackerHit ph = (*it)->first->edm4hepHit;
		    trk->addHit(ph) ;
		    ++nHit ;
		    //std::cout <<  "   hit  added  "<<  (*it)->first->edm4hepHit   << std::endl ;
		  }

		  trk->initialise( MarlinTrk::IMarlinTrack::backward ) ;
		}

		int code = trk->fit(  maxChi2  ) ;
                // for (auto hit : hitsInFit) std::cout << hit.first << std::endl;

		if( code != MarlinTrk::IMarlinTrack::success ){

			std::cout << "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IMarlinTrkFitter :  problem fitting track "
		          << " error code : " << MarlinTrk::errorCode( code )
		          << std::endl ;

		}


                std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hitsInFit ;
                trk->getHitsInFit( hitsInFit ) ;
		//----- if the fit did not fail but has a small number of hits used,
		//      we try again one more time with a larger max-chi2-increment
		if( isFirstFit  && ( 1.*hitsInFit.size()) / (1.*nHit )  <  0.2  ) {  // fixme: parameter

			isFirstFit = false ;

			maxChi2 =  2. * _maxChi2Increment  ;

			// streamlog_out( DEBUG4 ) << "  >>>>>>  IMarlinTrkFitter :  small number of hits used in fit " << hitsInFit.size() << "/" << nHit << " = "
			// << ( 1.*hitsInFit.size()) / (1.*nHit )  << " refit with larger max chi2 increment:  " << maxChi2 <<  std::endl ;
			delete trk ;

			goto start ;   // ;-)
		}
		//----------------------------------------------------------------------


		return trk;
	}


	//---------------------------------------------------------------------------------------------------------------------------

	edm4hep::Track PLCIOTrackConverter::operator() (CluTrack* c) {
	  
		static lcio::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;

		edm4hep::Track trk;

		trk.setType( lcio::ILDDetID::TPC ) ;

		double e = 0.0 ;
		int nHit = 0 ;
		for( CluTrack::iterator hi = c->begin(); hi != c->end() ; hi++) {

			// reset outliers (not used in fit)  bit
			//      IMPL::TrackerHitImpl* thi = dynamic_cast<IMPL::TrackerHitImpl*> (  (*hi)->first->edm4hepHit )  ;
			//      thi->setQualityBit( UTIL::ILDTrkHitQualityBit::USED_IN_FIT , 0 )  ;
		  
		        trk.addToTrackerHits(  (*hi)->first->edm4hepHit ) ;
			e += (*hi)->first->edm4hepHit.getEDep() ;
			nHit++ ;
		}

		MarlinTrk::IMarlinTrack* mtrk = MarTrkof(c);

		MarTrk_of_edm4hepTrack(trk) = mtrk ;

		trk.setDEdx( ( nHit ? e/nHit : -1. )  ) ;

		/* FIXME Mingrui
		   trk->subdetectorHitNumbers().resize( 2 * lcio::ILDDetID::ETD ) ;

		   trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 2 ] =  0 ;
		   trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 1 ] =  nHit ;
		   */

		RuntimeMap<edm4hep::ConstTrackerHit, int> DChi2_of_hit;

		if( mtrk != 0 && ! c->empty() ){


			std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hitsInFit ;
			mtrk->getHitsInFit( hitsInFit ) ;
                        // for (auto hit : hitsInFit) std::cout << hit.second << std::endl;
			// FIXME Mingrui
			// trk->subdetectorHitNumbers()[ 2*lcio::ILDDetID::TPC - 2 ] =  hitsInFit.size() ;

			if( ! hitsInFit.empty() ){

				// store the delta chi2 for the given hit
				for(unsigned i=0, N=hitsInFit.size() ; i<N ; ++i){

					DChi2_of_hit(hitsInFit[i].first) = hitsInFit[i].second ;

					//IMPL::TrackerHitImpl* thi = dynamic_cast<IMPL::TrackerHitImpl*> ( hitsInFit[i].first ) ;
					//thi->setQualityBit( UTIL::ILDTrkHitQualityBit::USED_IN_FIT , 1 )  ;
				}

				edm4hep::TrackState tsIP;
				edm4hep::TrackState tsFH;
				edm4hep::TrackState tsLH;
				edm4hep::TrackState tsCA;

				tsIP.location = lcio::TrackState::AtIP;
				tsFH.location = lcio::TrackState::AtFirstHit;
				tsLH.location = lcio::TrackState::AtLastHit;
				tsCA.location = lcio::TrackState::AtCalorimeter;

				double chi2 ;
				int ndf  ;
				int code ;

				// Hit* hf = c->front() ;
				// Hit* hb = c->back() ;
				// bool reverse_order =   ( std::abs( hf->first->pos.z() ) > std::abs( hb->first->pos.z()) + 3. ) ;
				// lcio::TrackerHit* fHit =  ( reverse_order ?  hb->first->lcioHit  :  hf->first->lcioHit ) ;
				// lcio::TrackerHit* lHit =  ( reverse_order ?  hf->first->lcioHit  :  hb->first->lcioHit ) ;

				edm4hep::ConstTrackerHit fHit = (hitsInFit.back().first);
				edm4hep::ConstTrackerHit lHit = (hitsInFit.front().first);

				//order of hits in fit is reversed wrt time  (we fit inwards)

				// ======= get TrackState at first hit  ========================

				code = mtrk->getTrackState( fHit, tsFH, chi2, ndf ) ;

				if( code != MarlinTrk::IMarlinTrack::success ){

					   std::cout << "  >>>>>>>>>>> PLCIOTrackConverter :  could not get TrackState at first Hit !!?? "
					   << " error code : " << MarlinTrk::errorCode( code )
					   << std::endl ;
				}

				// ======= get TrackState at last hit  ========================

#define use_fit_at_last_hit 0

#if use_fit_at_last_hit
				code = mtrk->getTrackState( lHit, tsLH, chi2, ndf ) ;
#else     // get the track state at the last hit by propagating from the last(first) constrained fit position (a la MarlinTrkUtils)
				edm4hep::ConstTrackerHit last_constrained_hit;
				code = mtrk->getTrackerHitAtPositiveNDF( last_constrained_hit );
				mtrk->smooth() ;
				if( code != MarlinTrk::IMarlinTrack::success ){
				  std::cout << "last_constrained_hit is unavaibale" << std::endl;
				}
				else{
// std::cout << "the hit" << std::endl;
// std::cout << lHit.getCellID0() << std::endl;
// std::cout << last_constrained_hit << std::endl;
				//code = mtrk->smooth() ;
// std::cout << "Smooth success" << std::endl;
				// gear::Vector3D last_hit_pos( lHit->getPosition()[0], lHit->getPosition()[1], lHit->getPosition()[2] );
// std::cout << "lHit = " << lHit << std::endl;

// std::cout << "----------------------------------------------------------" << std::endl;
// std::cout << "Position" << lHit.getPosition() << std::endl;
// std::cout << "----------------------------------------------------------" << std::endl;

				  edm4hep::Vector3d last_hit_pos( lHit.getPosition() );
// std::cout << "lHit = " << lHit << std::endl;
				  code = mtrk->propagate( last_hit_pos, last_constrained_hit, tsLH, chi2, ndf);
                                } 
// std::cout << "Propagate success" << std::endl;
#endif

// std::cout << "We are here" << std::endl; 
				if( code != MarlinTrk::IMarlinTrack::success ){

					std::cout << "  >>>>>>>>>>> PLCIOTrackConverter :  could not get TrackState at last Hit !!?? " << std::endl ;
				}

				// ======= get TrackState at calo face  ========================
				//
				encoder.reset() ;
				encoder[ lcio::ILDCellID0::subdet ] =  lcio::ILDDetID::ECAL ;
				encoder[ lcio::ILDCellID0::layer  ] =  0  ;
				encoder[ lcio::ILDCellID0::side   ] =  lcio::ILDDetID::barrel;
				int layerID  = encoder.lowWord() ;
				int sensorID = -1 ;

#if use_fit_at_last_hit
				code = mtrk->propagateToLayer( layerID , lHit, tsCA, chi2, ndf, sensorID, MarlinTrk::IMarlinTrack::modeClosest ) ;
#else     // get the track state at the calorimter from a propagating from the last(first) constrained fit position
				code = mtrk->propagateToLayer( layerID , last_constrained_hit, tsCA, chi2, ndf, sensorID, MarlinTrk::IMarlinTrack::modeClosest ) ;
#endif

// std::cout << "We are here 2" << std::endl; 
				if( code ==  MarlinTrk::IMarlinTrack::no_intersection ){

					encoder[ lcio::ILDCellID0::side   ] = ( lHit.getPosition()[2] > 0.  ?   lcio::ILDDetID::fwd  :  lcio::ILDDetID::bwd  ) ;
					layerID = encoder.lowWord() ;

#if use_fit_at_last_hit
					code = mtrk->propagateToLayer( layerID , lHit, tsCA, chi2, ndf, sensorID, MarlinTrk::IMarlinTrack::modeClosest ) ;
#else     // get the track state at the calorimter from a propagating from the last(first) constrained fit position
					code = mtrk->propagateToLayer( layerID , last_constrained_hit, tsCA, chi2, ndf, sensorID, MarlinTrk::IMarlinTrack::modeClosest ) ;
#endif

				}
				if ( code !=MarlinTrk::IMarlinTrack::success ) {

					// streamlog_out( DEBUG6 ) << "  >>>>>>>>>>> PLCIOTrackConverter :  could not get TrackState at calo face !!?? " << std::endl ;
				}

				// ======= get TrackState at IP ========================

// std::cout << "We are here 3" << std::endl; 
				// const gear::Vector3D ipv( 0.,0.,0. );
				const double tmp_ipv[] = {0, 0, 0};
				const edm4hep::Vector3d ipv(tmp_ipv);

				// fg: propagate is quite slow  and might not really be needed for the TPC

				code = ( UsePropagate ?   mtrk->propagate( ipv, fHit, tsIP, chi2, ndf ) :  mtrk->extrapolate( ipv, tsIP, chi2, ndf ) ) ;

				if( code != MarlinTrk::IMarlinTrack::success ){

					std::cout << "  >>>>>>>>>>> PLCIOTrackConverter :  could not extrapolate TrackState to IP !!?? " << std::endl ;
				}

				trk.addToTrackStates( tsIP ) ;
				trk.addToTrackStates( tsFH ) ;
				trk.addToTrackStates( tsLH ) ;
				trk.addToTrackStates( tsCA ) ;

// std::cout << "We are here 3" << std::endl; 
				double RMin = sqrt( tsFH.referencePoint[0] * tsFH.referencePoint[0]
						+ tsFH.referencePoint[1] * tsFH.referencePoint[1] ) ;

				trk.setRadiusOfInnermostHit( RMin  ) ;

				trk.setChi2( chi2 ) ;
				trk.setNdf( ndf ) ;

			} else {

				std::cout << "  >>>>>>>>>>> PLCIOTrackConverter::operator()  -  hitsInFitEmpty ! - nHits " << nHit << std::endl ;
			}

		} else {

			// this is not an error  for debug collections that just consist of track hits (no fit )
			// streamlog_out( ERROR ) << "  >>>>>>>>>>> LCIOTrackConverter::operator() (CluTrack* c)  :  "
			// 			     << " no MarlinTrk::IMarlinTrack* found for cluster !!?? " << std::endl ;
		}


		return trk;
	}


	}//namespace
