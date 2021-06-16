#include "ClupatraAlg.h"

#include "clupatra_new.h"

#include <time.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <memory>
#include <float.h>

//---- MarlinUtil

//---- LCIO ---
//#include "IMPL/LCCollectionVec.h"
//#include "EVENT/SimTrackerHit.h"
//#include "UTIL/Operators.h"
//#include "UTIL/LCTOOLS.h"
//#include "UTIL/CellIDDecoder.h"
#include "UTIL/ILDConf.h"
//#include "UTIL/LCIterator.h"

//-------gsl -----
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"


//---- GEAR ----
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/ZPlanarParameters.h"
#include "gear/ZPlanarLayerLayout.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"


#include "TrackSystemSvc/IMarlinTrack.h"
#include "TrackSystemSvc/IMarlinTrkSystem.h"
#include "TrackSystemSvc/MarlinTrkUtils.h"

#include "RuntimeMap.h"

#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
// #include "edm4hep/TrackerHitPlane.h"


using namespace edm4hep ;
using namespace MarlinTrk ;

using namespace clupatra_new ;

/*
   namespace edm4hep::TrackState {
   const int AtIP = 0;
   const int AtFirstHit = 1;
   const int AtLastHit = 2;
   const int AtCalorimeter = 3;
   };
   */

RuntimeMap<edm4hep::ConstTrack, clupatra_new::TrackInfoStruct*> TrackInfo_of_edm4hepTrack;
RuntimeMap<edm4hep::ConstTrack, MarlinTrk::IMarlinTrack*> MarTrk_of_edm4hepTrack;
RuntimeMap<MarlinTrk::IMarlinTrack*, clupatra_new::CluTrack*> CluTrk_of_MarTrack;
RuntimeMap<edm4hep::ConstTrackerHit, clupatra_new::Hit*> GHitof;
RuntimeMap<clupatra_new::CluTrack*, MarlinTrk::IMarlinTrack*> MarTrkof;

gear::GearMgr* gearMgr; 
#define WRITE_PICKED_DEBUG_TRACKS false

//----------------------------------------------------------------

template <class In, class Pred> In find_smallest(In first, In last, Pred p, double& d){

	In res =  last ;

	double min = DBL_MAX ;

	while( first != last ){

		double val = p( *first) ;

		if(  val < min ){

			res = first ;
			min = val ;
		}

		++first ;
	}

	d = min ;

	// FIXME: Mingrui ignore the debug

	return res ;
}
//----------------------------------------------------------------
struct Distance3D2{
	gear::Vector3D _pos ;
	Distance3D2( const gear::Vector3D& pos) : _pos( pos ) {}
	template <class T>
		double operator()( const T* t) {
			gear::Vector3D p( t->getPosition() ) ;
			return ( p - _pos ).r2() ;

		}
};

//----------------------------------------------------------------
struct StripDistance2{
	gear::Vector3D _pos ;
	StripDistance2( const gear::Vector3D& pos) : _pos( pos ) {}

	double operator()( const TrackerHit* t) {

/*
		gear::Vector3D p(t->getPosition()[0], t->getPosition()[1], t->getPosition()[2]) ;

		const TrackerHitPlane* h = (const TrackerHitPlane*) t ;

		gear::Vector3D v( 1. , h->getU()[1] ,  h->getU()[0] , gear::Vector3D::spherical ) ;

		double d = ( p - _pos ).dot(v)  ;

		return d*d ;
*/
return 1;
	}
};

//----------------------------------------------------------------
struct MeanAbsZOfTrack{
	double operator()( ConstTrack t){
		double z = 0 ;
		int hitCount = 0 ;
		/*
		   const TrackerHitVec& hV = t->getTrackerHits() ;
		   for(unsigned i=0, N = hV.size() ; i<N ; ++i){
		   z += hV[i]->getPosition()[2]  ;
		   ++hitCount ;
		   }
		   */
		for (auto iter = t.trackerHits_begin(); iter != t.trackerHits_end(); iter++) {
			z += iter->getPosition()[2];
			++hitCount;
		}
		return ( hitCount>0  ?  std::abs(  z )  / hitCount  : 0 ) ;
	}
};


DECLARE_COMPONENT( ClupatraAlg )

ClupatraAlg::ClupatraAlg(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc), _trksystem(0), _gearTPC(0) {

		// _description = "ClupatraProcessor : nearest neighbour clustering seeded pattern recognition" ;

		// Input Collections
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		declareProperty("TPCHitCollection", _TPCHitCollectionHandle, "handler of the tpc hit input collections");

		// Output Collections
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		declareProperty("ClupatraTracks", _ClupatraTrackCollectionHandle, "handler of the collection with final TPC tracks");
		declareProperty("ClupatraTrackSegments", _ClupatraTrackSegmentCollectionHandle, "handler of the collection that has the individual track segments");

	}

void ClupatraAlg::printParameters() {

  debug() << "Print parameters:" << endmsg;
  debug() << _distCut << endmsg;
  debug() << _cosAlphaCut << endmsg;
  debug() << _nLoop << endmsg;
  debug() << _minCluSize << endmsg;
  debug() << _duplicatePadRowFraction << endmsg; 
  debug() << _dChi2Max << endmsg;
  debug() << _chi2Cut << endmsg;
  debug() << _maxStep << endmsg;
  debug() << _pickUpSiHits << endmsg;
  debug() << _createDebugCollections << endmsg;
  debug() << _padRowRange << endmsg;
  debug() <<  _nZBins << endmsg;
  debug() << _minLayerFractionWithMultiplicity << endmsg;
  debug() << _minLayerNumberWithMultiplicity << endmsg;
  debug() << _trackStartsInnerDist << endmsg;
  debug() << _trackEndsOuterCentralDist << endmsg;
  debug() << _trackEndsOuterForwardDist << endmsg;
  debug() <<  _trackIsCurlerOmega << endmsg;

}


StatusCode ClupatraAlg::initialize() {

	// Usually a good idea to
        // don't need, since Gaudi Algorithm will print all Property  
	//printParameters() ;

	auto _trackSystemSvc = service<ITrackSystemSvc>("TrackSystemSvc");
	if ( !_trackSystemSvc ) {
		error() << "Fail to find TrackSystemSvc ..." << endmsg;
	}

	_trksystem = _trackSystemSvc->getTrackSystem(this);
	if(!_trksystem){
	  error() << "Cannot initialize MarlinTrkSystem of Type: KalTest" <<endmsg;
	  return StatusCode::FAILURE;
	}

	_trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
	_trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
	_trksystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
	_trksystem->init() ;

	_nRun = 0 ;
	_nEvt = 0 ;
	// FIXME: fucd
        //tree = new TTree("Tuple", "Particle Tree");
        //tree->Branch("omega", &omega, "omega/D");
        //tree->Branch("totalCandidates", &totalCandidates, "totalCandidates/I");
        //tree->Branch("eventNumber", &_nEvt, "eventNumber/I");

	return GaudiAlgorithm::initialize();

}


StatusCode ClupatraAlg::execute() {


	debug() << "Clupatra Algorithm started" << endmsg;

	//  clock_t start =  clock() ;
	Timer timer ;
	unsigned t_init       = timer.registerTimer(" initialization      " ) ;
	unsigned t_seedtracks = timer.registerTimer(" extend seed tracks  " ) ;
	unsigned t_recluster  = timer.registerTimer(" recluster leftovers " ) ;
	unsigned t_split      = timer.registerTimer(" split clusters      " ) ;
	unsigned t_finalfit   = timer.registerTimer(" final refit         " ) ;
	unsigned t_merge      = timer.registerTimer(" merge segments      " ) ;
	unsigned t_pickup     = timer.registerTimer(" pick up Si hits     " ) ;

	timer.start() ;

	// the clupa wrapper hits that hold pointers to PLCIO hits plus some additional parameters
	// create them in a vector for convenient memeory mgmt
	std::vector<ClupaHit> clupaHits ;

	// on top of the clupahits we need the tiny wrappers for clustering - they are created on the heap
	// and we put them in a vector of pointers that takes ownership (i.e. deletes them at the end)
	HitVec nncluHits ;
	nncluHits.setOwner( true ) ;


	// this is the final list of cluster tracks
	Clusterer::cluster_list cluList ;
	cluList.setOwner() ;


	PLCIOTrackConverter converter ;

	// Gear can be put as global
	auto _gear = service<IGearSvc>("GearSvc");
	gearMgr = _gear->getGearMgr();

	_gearTPC = &(gearMgr->getTPCParameters());
	_bfield = gearMgr->getBField().at( gear::Vector3D(0.,0.0,0.) ).z() ;

	// Support for more than one module
	// The ternary operator is used to make the trick with the static variable which
	// is supposed to be calculated only once, also for performance reason 
	static const unsigned int maxTPCLayers = (gearMgr->getDetectorName() == "LPTPC" ) ?
		_gearTPC->getModule(0).getNRows() + _gearTPC->getModule(2).getNRows() + _gearTPC->getModule(5).getNRows() :  // LCTPC
		_gearTPC->getModule(0).getNRows(); // ILD

	double driftLength = _gearTPC->getMaxDriftLength() ;
	ZIndex zIndex( -driftLength , driftLength , _nZBins  ) ;


	const edm4hep::TrackerHitCollection* col = nullptr;
        debug() << "col" << endmsg;

	try{   col = _TPCHitCollectionHandle.get();

	} catch(...) {
		// FIXME: Mingrui fix the output message
		// streamlog_out( WARNING ) <<  " input collection not in event : " << _colName << "   - nothing to do  !!! " << std::endl ;
	}

	//===============================================================================================
	//   create clupa and clustering hits for every edm4hep hit
	//===============================================================================================
	if (col == nullptr) {
	  return StatusCode::SUCCESS;
        }


	int nHit = col->size() ;

	clupaHits.resize( nHit ) ;       // creates clupa hits (w/ default c'tor)
	nncluHits.reserve( nHit ) ;

	debug()  << "  create clupatra TPC hits, n = " << nHit << endmsg ;

	// Mingrui!: Copy the items in col to col_copy to avoid losing the pointer.
	for(int i=0 ; i < nHit ; ++i ) {

		//------
		ConstTrackerHit th(col->at(i));
		//debug() << i << " " << th->getCellID() << endmsg;
		if ( fabs(th.getPosition()[2]) > driftLength ) continue;

		ClupaHit* ch  = & clupaHits[i] ;

		Hit* gh =  new Hit( ch ) ;

		nncluHits.push_back( gh ) ;

		//-------
		// FIXME: Here we should have a resolution
		GHitof(th) = gh ;  // assign the clupa hit to the LCIO hit for memory mgmt

		ch->edm4hepHit = th ;

		ch->pos = gear::Vector3D(th.getPosition()[0], th.getPosition()[1], th.getPosition()[2]) ;

		//  int padIndex = padLayout.getNearestPad( ch->pos.rho() , ch->pos.phi() ) ;
		//    ch->layer = padLayout.getRowNumber( padIndex ) ;
		// May cause a problem
		ch->layer = getLayer( th );

		//debug() << "ch->layer = idDec( th )[ ILDCellID0::layer ] = " <<  ch->layer << " - CellID0 " << th.getCellID() << endmsg;

		ch->zIndex = zIndex( &th ) ;

		//ch->phiIndex = ....

	}
        debug() << "Get hits for clustering" << endmsg;

	//---------------------------------------------------------------------------------------------------------

	std::sort( nncluHits.begin(), nncluHits.end() , ZSort() ) ;

	//---------------------------------------------------------------------------------------------------------

	HitListVector hitsInLayer( maxTPCLayers ) ;
	addToHitListVector(  nncluHits.begin(), nncluHits.end() , hitsInLayer  ) ;

	debug() << "  added  " <<  nncluHits.size()  << "  tp hitsInLayer - > size " <<  hitsInLayer.size() << endmsg;

	//---------------------------------------------------------------------------------------------------------

	//===============================================================================================
	//   create output collections  ( some optional )
	//===============================================================================================

	// FIXME: Mingrui how to create the correct output collection?

	// const bool writeSeedCluster        = _createDebugCollections ;
	// const bool writeCluTrackSegments   = _createDebugCollections ;
	// const bool writeLeftoverClusters   = _createDebugCollections ;
	// const bool writeQualityTracks      = _createDebugCollections ;
	// const bool writeDebugTracks      =  WRITE_PICKED_DEBUG_TRACKS ;

	static const bool copyTrackSegments = false ;

	/*
	 * FIXME Mingrui remove these collections for now
	 LCCollectionVec* seedCol =  ( writeSeedCluster        ?  newTrkCol( "ClupatraSeedCluster"          , evt )  :   0   )  ;
	 LCCollectionVec* cluCol  =  ( writeCluTrackSegments   ?  newTrkCol( "ClupatraInitialTrackSegments" , evt )  :   0   )  ;
	 LCCollectionVec* locCol  =  ( writeCluTrackSegments   ?  newTrkCol( "ClupatraLeftoverClusters"     , evt )  :   0   )  ;

	 LCCollectionVec* incSegCol  = ( writeCluTrackSegments   ?  newTrkCol( "ClupatraIncompleteSegments"   , evt , true )  :   0   )  ;
	 LCCollectionVec* curSegCol  = ( writeCluTrackSegments   ?  newTrkCol( "ClupatraCurlerSegments"       , evt , true )  :   0   )  ;
	 LCCollectionVec* finSegCol  = ( writeCluTrackSegments   ?  newTrkCol( "ClupatraFinalTrackSegments"   , evt , true )  :   0   )  ;

	//LCCollectionVec* goodCol  = ( writeQualityTracks ?  newTrkCol( "ClupatraGoodQualityTracks" , evt ,true )  :   0   )  ;
	//LCCollectionVec* fairCol  = ( writeQualityTracks ?  newTrkCol( "ClupatraFairQualityTracks" , evt ,true )  :   0   )  ;
	LCCollectionVec* poorCol  = ( writeQualityTracks ?  newTrkCol( "ClupatraPoorQualityTracks" , evt , true )  :   0   )  ;

	LCCollectionVec* debugCol=  ( writeDebugTracks ?  newTrkCol( "ClupatraDebugTracks" , evt , false )  :   0   )  ;
	if( WRITE_PICKED_DEBUG_TRACKS )
	DebugTracks::setCol( debugCol , this ) ;

	LCCollectionVec* outerCol  = ( _createDebugCollections ?  newTrkCol( "ClupatraOuterSegments" , evt ,true )  :   0   )  ;
	LCCollectionVec* innerCol  = ( _createDebugCollections ?  newTrkCol( "ClupatraInnerSegments" , evt ,true )  :   0   )  ;
	LCCollectionVec* middleCol = ( _createDebugCollections ?  newTrkCol( "ClupatraMiddleSegments" , evt ,true )  :   0   )  ;
	*/



	TrackCollection* outCol =  _ClupatraTrackCollectionHandle.createAndPut();
	TrackCollection* tsCol  =  _ClupatraTrackSegmentCollectionHandle.createAndPut();
        std::vector<ClupaPlcioTrack*> outCol_tmp;
        std::vector<ClupaPlcioTrack*> tsCol_tmp;
	
	//---------------------------------------------------------------------------------------------------------

	timer.time(t_init ) ;

	//===============================================================================================
	// first main step of clupatra:
	//   * cluster in pad row range - starting from the outside - to find clean cluster segments
	//   * extend the track segments with best matching hits, based on extrapolation to next layer(s)
	//   * add the hits and apply a Kalman filter step ( track segement is always best estimate )
	//   * repeat in backward direction ( after smoothing back, to get a reasonable track
	//     state for extrapolating backwards )
	//===============================================================================================

	Clusterer nncl ;

	int outerRow = 0 ;

	nnclu::PtrVector<MarlinTrk::IMarlinTrack> seedTrks ;
	seedTrks.setOwner() ; // memory mgmt - will delete MarlinTrks at the end

	IMarlinTrkFitter fitter( _trksystem ) ;


	debug() << "===============================================================================================" << endmsg;
	debug() << "   first step of Clupatra algorithm: find seeds with NN-clustering  in " <<  _nLoop << " loops - max dist = " << _distCut << endmsg;
	debug() << "==============================================================================================="  << endmsg;

	// ---- introduce a loop over increasing distance cuts for finding the tracks seeds
	//      -> should fix (some of) the problems seen @ 3 TeV with extremely boosted jets
	//
	double dcut =  _distCut / _nLoop ;
	for(int nloop=1 ; nloop <= _nLoop ; ++nloop) {

		HitDistance dist( nloop * dcut , _cosAlphaCut ) ;

		outerRow = maxTPCLayers - 1 ;

		while( outerRow >= _minCluSize ) {

			HitVec hits ;
			hits.reserve( nHit ) ;

			// add all hits in pad row range to hits
			for(int iRow = outerRow ; iRow > ( outerRow - _padRowRange) ; --iRow ) {

				if( iRow > -1 ) {

					// debug() << "  copy " <<  hitsInLayer[ iRow ].size() << " hits for row " << iRow << endmsg ;

					std::copy( hitsInLayer[ iRow ].begin() , hitsInLayer[ iRow ].end() , std::back_inserter( hits )  ) ;
				}
			}

			//-----  cluster in given pad row range  -----------------------------
			Clusterer::cluster_list sclu ;
			sclu.setOwner() ;

			// FIXME Mingrui
			debug() << "   call cluster_sorted with " <<  hits.size() << " hits " << endmsg;

			nncl.cluster_sorted( hits.begin(), hits.end() , std::back_inserter( sclu ), dist , _minCluSize ) ;

			const static int merge_seeds = true ;

			if( merge_seeds ) { //-----------------------------------------------------------------------

				// sometimes we have split seed clusters as one link is just above the cut
				// -> recluster in all hits of small clusters with 1.2 * cut
				float _smallClusterPadRowFraction = 0.9  ;
				float _cutIncrease = 1.2 ;
				// fixme: could make parameters ....

				HitVec seedhits ;
				Clusterer::cluster_list smallclu ;
				smallclu.setOwner() ;
				split_list( sclu, std::back_inserter(smallclu),  ClusterSize(  int( _padRowRange * _smallClusterPadRowFraction) ) ) ;
				for( Clusterer::cluster_list::iterator sci=smallclu.begin(), end= smallclu.end() ; sci!=end; ++sci ){
					for( Clusterer::cluster_type::iterator ci=(*sci)->begin(), end= (*sci)->end() ; ci!=end;++ci ){
						seedhits.push_back( *ci ) ;
					}
				}
				// free hits from bad clusters
				std::for_each( smallclu.begin(), smallclu.end(), std::mem_fun( &CluTrack::freeElements ) ) ;

				HitDistance distLarge( nloop * dcut * _cutIncrease ) ;

				nncl.cluster_sorted( seedhits.begin(), seedhits.end() , std::back_inserter( sclu ), distLarge , _minCluSize ) ;

			} //------------------------------------------------------------------------------------------

			// FIXME Mingrui
			debug()  << "     found " <<  sclu.size() << "  clusters " << endmsg;

			// try to split up clusters according to multiplicity
			int layerWithMultiplicity = _padRowRange - 2  ; // fixme: make parameter
			split_multiplicity( sclu , layerWithMultiplicity , 10 ) ;


			// remove clusters whith too many duplicate hits per pad row
			Clusterer::cluster_list bclu ;    // bad clusters
			bclu.setOwner() ;
			split_list( sclu, std::back_inserter(bclu),  DuplicatePadRows( maxTPCLayers, _duplicatePadRowFraction  ) ) ;
			// free hits from bad clusters
			std::for_each( bclu.begin(), bclu.end(), std::mem_fun( &CluTrack::freeElements ) ) ;


			// ---- now we also need to remove the hits from good cluster seeds from the hitsInLayers:
			for( Clusterer::cluster_list::iterator sci=sclu.begin(), end= sclu.end() ; sci!=end; ++sci ){
				for( Clusterer::cluster_type::iterator ci=(*sci)->begin(), end= (*sci)->end() ; ci!=end;++ci ){

					// this is not cheap ...
					hitsInLayer[ (*ci)->first->layer ].remove( *ci )  ;
				}
			}

			// now we have 'clean' seed clusters
			// Write debug collection with seed clusters:
			// convert the clusters into tracks and write the resulting tracks into the existing debug collection seedCol.
			// The conversion is performed by the STL transform() function, the insertion to the end of the
			// debug track collection is done by creating an STL back_inserter iterator on the LCCollectionVector seedCol
			/* Debug FIXME Mingrui
			   if( writeSeedCluster ) {
			   std::transform( sclu.begin(), sclu.end(), std::back_inserter( *seedCol ) , converter ) ;
			   }
			   */

			//      std::transform( sclu.begin(), sclu.end(), std::back_inserter( seedTrks) , fitter ) ;
			// reduce memory footprint: deal with one KalTest track at a time and delete it, when done

			for( Clusterer::cluster_list::iterator icv = sclu.begin(), end =sclu.end()  ; icv != end ; ++ icv ) {
				// icv = *cluster_list type

				int nHitsAdded = 0 ;

				debug() <<  " call fitter for seed cluster with " << (*icv)->size() << " hits " << endmsg;
				int counter = 0;
				for( Clusterer::cluster_type::iterator ci=(*icv)->begin(), end= (*icv)->end() ; ci!=end; ++ci ) {
				  debug() << counter++ << " " <<  (*ci)->first->edm4hepHit << " \nlayer " << (*ci)->first->layer << endmsg;
				}


				MarlinTrk::IMarlinTrack* mTrk = fitter( *icv ) ;
				debug() << "before add hits and filter" << endmsg;
                // std::vector<std::pair<edm4hep::ConstTrackerHit, double> > hitsInFit ;
                // mTrk->getHitsInFit( hitsInFit ) ;
                // for (auto hit : hitsInFit) std::cout << hit.first << std::endl;

				nHitsAdded += addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex ) ;

				static const bool backward = true ;
				nHitsAdded += addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ;
				// in order to use smooth for backward extrapolation call with   _trksystem  - does not work well...
				// nHitsAdded += addHitsAndFilter( *icv , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward , _trksystem ) ;


				// drop seed clusters with no hits added - but not in the very forward region...
				debug() << "Goes here" << endmsg;
				if( nHitsAdded < 1  &&  outerRow >   2*_padRowRange  ){  //FIXME: make parameter ?

					ConstTrack edm4hepTrk( converter( *icv ) ) ;
					// debug() << "Goes goes here" << endmsg;

					// FIXME Mingrui
					debug() << "=============  poor seed cluster - no hits added - started from row " <<  outerRow << "\n"
						<< edm4hepTrk << endmsg;


					for( Clusterer::cluster_type::iterator ci=(*icv)->begin(), end= (*icv)->end() ; ci!=end; ++ci ) {
						hitsInLayer[ (*ci)->first->layer ].push_back( *ci )   ;
					}
					(*icv)->freeElements() ;
					(*icv)->clear() ;
				}

				/*
				   if( writeCluTrackSegments )  //  ---- store track segments from the first main step  -----
				   cluCol->addElement(  converter( *icv ) );
				   Debug FIXME Mingrui
				   */
				// reset the pointer to the KalTest track - as we are done with this track
				MarTrkof(*icv) = 0;

				delete mTrk ;
			}

			// merge the good clusters to final list
			cluList.merge( sclu ) ;

			outerRow -= _padRowRange ;

		} //while outerRow > padRowRange

	} // nloop


	timer.time( t_seedtracks ) ;

	timer.time( t_recluster ) ;

	//===============================================================================================
	//  do a global reclustering in leftover hits
	//===============================================================================================
	debug() << "do global reclustering" << endmsg;
	static const int do_global_reclustering = true ;
	if( do_global_reclustering ) {

		outerRow = maxTPCLayers - 1 ;

		int padRangeRecluster = 50 ; // FIXME: make parameter
		// define an inner cylinder where we exclude hits from re-clustering:
		double zMaxInnerHits   = driftLength * .67 ;   // FIXME: make parameter
		double rhoMaxInnerHits = _gearTPC->getPlaneExtent()[0] + (  _gearTPC->getPlaneExtent()[1] - _gearTPC->getPlaneExtent()[0] ) * .67 ;// FIXME: make parameter


		while( outerRow > 0 ) {


			Clusterer::cluster_list loclu ; // leftover clusters
			loclu.setOwner() ;

			HitVec hits ;

			int  minRow = ( ( outerRow - padRangeRecluster ) > -1 ?  ( outerRow - padRangeRecluster ) : -1 ) ;

			// add all hits in pad row range to hits
			for(int iRow = outerRow ; iRow > minRow ; --iRow ) {


				for( HitList::iterator hlIt=hitsInLayer[ iRow ].begin() , end = hitsInLayer[ iRow ].end() ; hlIt != end ; ++hlIt ) {

					if( std::abs( (*hlIt)->first->pos.z() ) > zMaxInnerHits  ||  (*hlIt)->first->pos.rho() >  rhoMaxInnerHits ) {
						hits.push_back( *hlIt ) ;
					}
				}
			}


			HitDistance distSmall( _distCut ) ;
			nncl.cluster( hits.begin(), hits.end() , std::back_inserter( loclu ),  distSmall , _minCluSize ) ;

			// Write debug collection using STL transform() function on the clusters
			/* Debug FIXME Mingrui
			   if( writeLeftoverClusters )
			   std::transform( loclu.begin(), loclu.end(), std::back_inserter( *locCol ) , converter ) ;
			   */


			// timer.time( t_recluster ) ;

			//===============================================================================================
			//  now we split the clusters based on their hit multiplicities
			//===============================================================================================


			//    _dChi2Max = 5. * _dChi2Max ; //FIXME !!!!!!!!!

			for( Clusterer::cluster_list::iterator it= loclu.begin(), end= loclu.end() ; it != end ; ++it ){

				CluTrack* clu = *it ;

				// FIXME: Mingrui

				std::vector<int> mult(8) ;
				// get hit multiplicities up to 6 ( 7 means 7 or higher )
				getHitMultiplicities( clu , mult ) ;

				// FIXME Mingrui
				// streamlog_out(  DEBUG3 ) << " **** left over cluster with hit multiplicities: \n" ;
				for( unsigned i=0,n=mult.size() ; i<n ; ++i) {
					// FIXME Mingrui
					// streamlog_out(  DEBUG3 ) << "     m["<<i<<"] = " <<  mult[i] << "\n"  ;
				}


				if( float( mult[5]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[5] >  _minLayerNumberWithMultiplicity ) {

					Clusterer::cluster_list reclu ; // reclustered leftover clusters
					reclu.setOwner() ;

					create_n_clusters( *clu , reclu , 5 ) ;

					std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;

					for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){

						// FIXME Mingrui
						// streamlog_out( DEBUG5 ) << " extending mult-5 clustre  of length " << (*ir)->size() << std::endl ;

						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ;
						static const bool backward = true ;
						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ;
					}

					cluList.merge( reclu ) ;
				}

				else if( float( mult[4]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[4] >  _minLayerNumberWithMultiplicity ) {

					Clusterer::cluster_list reclu ; // reclustered leftover clusters
					reclu.setOwner() ;

					create_n_clusters( *clu , reclu , 4 ) ;

					std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;

					for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){


						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ;
						static const bool backward = true ;
						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ;
					}

					cluList.merge( reclu ) ;
				}

				else if( float( mult[3]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[3] >  _minLayerNumberWithMultiplicity ) {

					Clusterer::cluster_list reclu ; // reclustered leftover clusters
					reclu.setOwner() ;

					create_three_clusters( *clu , reclu ) ;

					std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;

					for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){

						// FIXME Mingrui
						// streamlog_out( DEBUG5 ) << " extending triplet clustre  of length " << (*ir)->size() << std::endl ;

						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ;
						static const bool backward = true ;
						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ;
					}

					cluList.merge( reclu ) ;
				}

				else if( float( mult[2]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[2] >  _minLayerNumberWithMultiplicity ) {

					Clusterer::cluster_list reclu ; // reclustered leftover clusters
					reclu.setOwner() ;

					create_two_clusters( *clu , reclu ) ;

					std::transform( reclu.begin(), reclu.end(), std::back_inserter( seedTrks) , fitter ) ;

					for( Clusterer::cluster_list::iterator ir= reclu.begin(), end= reclu.end() ; ir != end ; ++ir ){

						// FIXME Mingrui
						// streamlog_out( DEBUG5 ) << " extending doublet clustre  of length " << (*ir)->size() << std::endl ;

						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ;
						static const bool backward = true ;
						addHitsAndFilter( *ir , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ;
					}

					cluList.merge( reclu ) ;

				}
				else if( float( mult[1]) / mult[0]  >= _minLayerFractionWithMultiplicity &&  mult[1] >  _minLayerNumberWithMultiplicity ) {


					seedTrks.push_back( fitter( *it )  );

					addHitsAndFilter( *it , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex) ;
					static const bool backward = true ;
					addHitsAndFilter( *it , hitsInLayer , _dChi2Max, _chi2Cut , _maxStep , zIndex, backward ) ;

					cluList.push_back( *it ) ;

					it = loclu.erase( it ) ;
					--it ; // erase returns iterator to next element

				} else {

					//  discard cluster and free hits
					clu->freeElements() ;
				}

			}


			outerRow -=  padRangeRecluster ;

		}
	}

	//=======================================================================================================================
	//  try again to gobble up hits at the ends ....   - does not work right now, as there are no fits  for the clusters....
	//=======================================================================================================================


	timer.time( t_split ) ;

	//===============================================================================================
	//  now refit the tracks
	//===============================================================================================

	// FIXME Mingrui
	debug()  << " ===========    refitting final " << cluList.size() << " track segments  "   << endmsg ;

	//---- refit cluster tracks individually to save memory ( KalTest tracks have ~1MByte each)

	IMarlinTrkFitter fit(_trksystem,  _dChi2Max) ; // fixme: do we need a different chi2 max here ????

	for( Clusterer::cluster_list::iterator icv = cluList.begin() , end = cluList.end() ; icv != end ; ++ icv ) {

		if( (*icv)->empty() )
			continue ;

		MarlinTrk::IMarlinTrack* trk = fit( *icv ) ;
		trk->smooth() ;
		edm4hep::Track edm4hepTrk = converter( *icv ) ;
		tsCol_tmp.push_back( new ClupaPlcioTrack(edm4hepTrk) ) ;
		MarTrk_of_edm4hepTrack(edm4hepTrk) = 0 ;
		delete trk ;
	}

	timer.time( t_finalfit) ;

	//===============================================================================================
	//   optionally create collections of used and unused TPC hits
	//===============================================================================================


	//===============================================================================================
	//  compute some track parameters for possible merging
	//===============================================================================================

	typedef nnclu::NNClusterer<ClupaPlcioTrack> TrackClusterer ;
	TrackClusterer nntrkclu ;
	MakePLCIOElement<ClupaPlcioTrack> trkMakeElement ;

	for( int i=0,N=tsCol_tmp.size() ;  i<N ; ++i ) {
	  edm4hep::ConstTrack track = tsCol_tmp.at(i)->edm4hepTrack;
	  computeTrackInfo(track) ;
	}

	//===============================================================================================
	//  merge split segements
	//===============================================================================================


	static const int merge_split_segments = true ;

	if( merge_split_segments ) {

		for(unsigned l=0 ; l < 2 ; ++l ) { // do this twice ....

			int nMax  =  tsCol_tmp.size()   ;

			TrackClusterer::element_vector incSegVec ;
			incSegVec.setOwner() ;
			incSegVec.reserve( nMax  ) ;
			TrackClusterer::cluster_vector incSegCluVec ;
			incSegCluVec.setOwner() ;

			for( int i=0,N=tsCol_tmp.size() ;  i<N ; ++i ){

			        edm4hep::ConstTrack trk = tsCol_tmp.at(i)->edm4hepTrack;

				const TrackInfoStruct* ti = TrackInfo_of_edm4hepTrack(trk);

				bool isIncompleteSegment =   !ti->isCurler  && ( !ti->startsInner || ( !ti->isCentral && !ti->isForward )  ) ;

				std::bitset<32> type = trk.getType() ;

				if( isIncompleteSegment  && ! type[ lcio::ILDTrackTypeBit::SEGMENT ]){

					incSegVec.push_back(  trkMakeElement( tsCol_tmp[i] )  ) ;

					// FIXME Mingrui Debug part
					/*
					   if( writeCluTrackSegments )  incSegCol->addElement( trk ) ;
					   */
				}
			}


			TrackSegmentMerger trkMerge( _dChi2Max , _trksystem ,  _bfield  ) ;

			nntrkclu.cluster( incSegVec.begin() , incSegVec.end() , std::back_inserter( incSegCluVec ), trkMerge , 2  ) ;

			// FIXME: Mingrui
			// streamlog_out( DEBUG4 ) << " ===== merged track segments - # cluster: " << incSegCluVec.size()
			//  << " from " << incSegVec.size() << " incomplete track segments "    << "  ============================== " << std::endl ;

			for(  TrackClusterer::cluster_vector::iterator it= incSegCluVec.begin() ; it != incSegCluVec.end() ; ++it) {

				// FIXME: Mingrui
				// streamlog_out( DEBUG4 ) <<  edm4hep::header<edm4hep::Track>() << std::endl ;

				TrackClusterer::cluster_type*  incSegClu = *it ;

				std::vector<edm4hep::ConstTrack> mergedTrk ;

				// vector to collect hits from segments
				//      std::vector< TrackerHit* >  hits ;
				// hits.reserve( 1024 ) ;
				// IMPL::TrackImpl* track = new  IMPL::TrackImpl ;
				// tsCol->addElement( track ) ;

				CluTrack hits ;

				for( TrackClusterer::cluster_type::iterator itC = incSegClu->begin() ; itC != incSegClu->end() ; ++ itC ){

					//streamlog_out( DEBUG3 ) << lcshort(  (*itC)->first ) << std::endl ;

					edm4hep::ConstTrack trk = (*itC)->first->edm4hepTrack;

					mergedTrk.push_back( trk ) ;

					//	std::copy( trk->getTrackerHits().begin() , trk->getTrackerHits().end() , std::back_inserter( hits ) ) ;

					/*
					   for( edm4hep::TrackerHitVec::const_iterator it = trk->getTrackerHits().begin() , END =  trk->getTrackerHits().end() ; it != END ; ++ it ){
					   hits.addElement( GHitof(*it) )  ;
					   }
					   */
					for (auto it = trk.trackerHits_begin(); it != trk.trackerHits_end(); it++) {
						ConstTrackerHit hit = *it;
						hits.addElement( GHitof(hit) );
					}

					// flag the segments so they can be ignored for final list
					// FIXME I need to set type
					// trk->setType( lcio::ILDTrackTypeBit::SEGMENT ) ;

					// add old segments to new track
					//	track->addTrack( trk ) ;
				}

				// MarlinTrk::IMarlinTrack* mTrk = _trksystem->createTrack();
				// EVENT::FloatVec icov( 15 ) ;
				// icov[ 0] = 1e2 ;
				// icov[ 2] = 1e2 ;
				// icov[ 5] = 1e2 ;
				// icov[ 9] = 1e2 ;
				// icov[14] = 1e2 ;
				// int result = createFinalisedLCIOTrack( mTrk, hits, track, !MarlinTrk::IMarlinTrack::backward, icov, _bfield,  _dChi2Max ) ;
				// //int result = createFinalisedLCIOTrack( mTrk, hits, track, ! MarlinTrk::IMarlinTrack::backward, icov, _bfield,  _dChi2Max ) ;
				// // ???

				MarlinTrk::IMarlinTrack* mTrk = fit( &hits ) ;
				mTrk->smooth() ;
				edm4hep::Track track = converter( &hits ) ;
				tsCol_tmp.push_back( new ClupaPlcioTrack(track) ) ;
				MarTrk_of_edm4hepTrack(track) = 0 ;
				delete mTrk ;
				computeTrackInfo( track ) ;

				// FIXME: Mingrui
				// streamlog_out( DEBUG4 ) << "   ******  created new track : " << " : " << lcshort( (Track*) track )  << std::endl ;

			}
		}// loop over l
	}
	//===============================================================================================
	//  merge curler segments
	//===============================================================================================


	static const int merge_curler_segments = true ;

	if( merge_curler_segments ) {

		int nMax  =  tsCol_tmp.size()   ;

		TrackClusterer::element_vector curSegVec ;
		curSegVec.setOwner() ;
		curSegVec.reserve( nMax  ) ;
		TrackClusterer::cluster_vector curSegCluVec ;
		curSegCluVec.setOwner() ;

		//for(std::vector<ClupaPlcioTrack*>::iterator it=tsCol_tmp.begin();it!=tsCol_tmp.end();it++){
		for( int i=tsCol_tmp.size()-1 ;  i>=0 ; --i ){

		        edm4hep::Track trk = tsCol_tmp.at(i)->edm4hepTrack;

			std::bitset<32> type = trk.getType() ;

			if( type[ lcio::ILDTrackTypeBit::SEGMENT ] )
				continue ;   // ignore previously merged track segments

			const TrackInfoStruct* ti = TrackInfo_of_edm4hepTrack(trk) ;

			bool isCompleteTrack =   ti && !ti->isCurler  && ( ti->startsInner &&  (  ti->isCentral || ti->isForward ) );

			if( !isCompleteTrack ){

				curSegVec.push_back(  trkMakeElement( tsCol_tmp[i] )  ) ;

				// FIXME
				// if( writeCluTrackSegments )  curSegCol->push_back( trk ) ;

			} else {   // ... is not a curler ->  add a copy to the final tracks collection


			        if( copyTrackSegments) {
				  
				        outCol_tmp.push_back( new ClupaPlcioTrack(trk)  ) ;

				}else{

				        outCol_tmp.push_back( new ClupaPlcioTrack(trk) ) ;
					debug() << "now outCol_tmp size = " << outCol_tmp.size() << endmsg;
					// FIXME !!! Possibly different with the LCIO version
					// tsCol-> Do not remove it at the Beginning
					// Very very important
					// tsCol->removeElementAt( i ) ;
					// fucd
					//delete tsCol_tmp.at(i);
					tsCol_tmp.erase(tsCol_tmp.begin()+i);
				}

				// FIXME Mingrui remove
				// if( writeCluTrackSegments )  finSegCol->push_back( trk ) ;
			}
		}

		//======================================================================================================


		TrackCircleDistance trkMerge( 0.1 ) ;

		nntrkclu.cluster( curSegVec.begin() , curSegVec.end() , std::back_inserter( curSegCluVec ), trkMerge , 2  ) ;

		debug() << " ===== merged tracks - # cluster: " << curSegCluVec.size()
			<< " from " << tsCol_tmp.size() << " track segments "    << "  ============================== " << endmsg ;

		//for (auto t : tsCol_tmp)
		//debug() << t->edm4hepTrack << endmsg;

// debug() << "Finish this line" << endmsg;

		for(  TrackClusterer::cluster_vector::iterator it= curSegCluVec.begin() ; it != curSegCluVec.end() ; ++it) {

			// FIXME: Mingrui
			// streamlog_out( DEBUG4 ) <<  edm4hep::header<Track>() << std::endl ;

			TrackClusterer::cluster_type*  curSegClu = *it ;

			std::list<edm4hep::Track> mergedTrk ;

			for( TrackClusterer::cluster_type::iterator itC = curSegClu->begin() ; itC != curSegClu->end() ; ++ itC ){

			  //debug() << lcshort(  (*itC)->first ) << endmsg;
			  //debug() << getOmega((*itC)->first->edm4hepTrack) << endmsg;
			  mergedTrk.push_back( (*itC)->first->edm4hepTrack ) ;
			}


			mergedTrk.sort( TrackZSort() ) ;

			//================================================================================


			if( copyTrackSegments) {

				// ====== create a new LCIO track for the merged cluster ...
				// edm4hep::Track trk;
				edm4hep::Track trk; 

				trk.setType( lcio::ILDDetID::TPC ) ;

				// == and copy all the hits
				unsigned hitCount = 0 ;
				for( std::list<edm4hep::Track>::iterator itML = mergedTrk.begin() ; itML != mergedTrk.end() ; ++ itML ){

					for (auto itHit = (*itML).trackerHits_begin(); itHit != (*itML).trackerHits_end(); itHit++) {
						trk.addToTrackerHits(*itHit);
					}
					hitCount  += (*itML).trackerHits_size()  ;

					// add a pointer to the original track segment
					trk.addToTracks(*itML) ;
				}

				// take track states from first and last track :
				ConstTrack firstTrk = mergedTrk.front() ;
				ConstTrack lastTrk  = mergedTrk.back() ;

				/* !!!!!!!!!!!!!! critical important FIXME should wait zoujiaheng
				   edm4hep::TrackState ts;
				   ts = firstTrk->getTrackState( lcio::TrackState::AtIP  ) ;
				   if( ts ) trk->addTrackState( new TrackState( *ts )  ) ;

				   ts = firstTrk->getTrackState( lcio::TrackState::AtFirstHit  ) ;
				   if( ts ) 	trk->addTrackState( new TrackState( *ts )  ) ;

				   ts = lastTrk->getTrackState( lcio::TrackState::AtLastHit  ) ;
				   if( ts ) trk->addTrackState( new TrackState( *ts )  ) ;

				   ts = lastTrk->getTrackState( lcio::TrackState::AtCalorimeter  ) ;
				   if( ts ) trk->addTrackState( new TrackState( *ts )  ) ;
				   */
				if (hasTrackStateAt(firstTrk, lcio::TrackState::AtIP )) trk.addToTrackStates(getTrackStateAt(firstTrk, lcio::TrackState::AtIP));
				if (hasTrackStateAt(firstTrk, lcio::TrackState::AtFirstHit )) trk.addToTrackStates(getTrackStateAt(firstTrk, lcio::TrackState::AtFirstHit));
				if (hasTrackStateAt(lastTrk, lcio::TrackState::AtLastHit )) trk.addToTrackStates(getTrackStateAt(lastTrk, lcio::TrackState::AtLastHit));
				if (hasTrackStateAt(lastTrk, lcio::TrackState::AtCalorimeter )) trk.addToTrackStates(getTrackStateAt(lastTrk, lcio::TrackState::AtCalorimeter));


				MarTrk_of_edm4hepTrack(trk) = MarTrk_of_edm4hepTrack(firstTrk) ;

				// FIXME: Mingrui Maybe this info not needed
				//int hitsInFit  =  firstTrk->getSubdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ] ;
				trk.setChi2(     firstTrk.getChi2()     ) ;
				trk.setNdf(      firstTrk.getNdf()      ) ;
				trk.setDEdx(     firstTrk.getDEdx()     ) ;
				trk.setDEdxError(firstTrk.getDEdxError()) ;

				/* FIXME: Mingrui Maybe this subdetectorHitNumbers is not needed?
				   trk->subdetectorHitNumbers().resize( 2 * lcio::ILDDetID::ETD ) ;
				   trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ] =  hitsInFit ;
				   trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ] =  hitCount ;
				   */

				double RMin = 0.0;
				if (hasTrackStateAt(trk, lcio::TrackState::AtFirstHit  ) ) {
					TrackState ts = getTrackStateAt(trk, lcio::TrackState::AtFirstHit  ) ;
					RMin = sqrt( ts.referencePoint[0] * ts.referencePoint[0]
							+ ts.referencePoint[1] * ts.referencePoint[1] );

				}
				/*
				   double RMin = ( ts ?
				   sqrt( ts->referencePoint[0] * ts->referencePoint[0]
				   + ts->referencePoint[1] * ts->referencePoint[1] )
				   :  0.0 ) ;
				   */
				trk.setRadiusOfInnermostHit( RMin  ) ;



				//debug() << "   create new merged track from bestTrack parameters - ptr to MarlinTrk : " << trk->ext<MarTrk>()
				//	<< "   with subdetector hit numbers  : " <<  trk->subdetectorHitNumbers()[0 ] << " , " <<  trk->subdetectorHitNumbers()[1]
				//	<< endmsg;


				outCol_tmp.push_back( new ClupaPlcioTrack(trk) )  ;

			} else { //==========================

				// we move the first segment to the final list and keep pointers to the other segments

			        std::list<edm4hep::Track>::iterator itML = mergedTrk.begin() ;

				edm4hep::Track trk = (edm4hep::Track) *itML++ ;

				for(  ; itML != mergedTrk.end() ; ++itML ){

					// add a pointer to the original track segment
					// FIXME I need to add sub track
				  trk.addToTracks(*itML) ;
				}

				outCol_tmp.push_back( new ClupaPlcioTrack(trk) ) ;
				debug() << "now outCol_tmp size = " << outCol_tmp.size() << endmsg;
				//remove from segment collection:
				for( int i=0,N=tsCol_tmp.size() ; i<N ; ++i) {
				  if( ( tsCol_tmp.at(i)->edm4hepTrack ) == trk ){
					        // FIXME Do not remove
						//tsCol_tmp.removeElementAt( i ) ;
					  //fucd
					  //delete tsCol_tmp[i];
					  tsCol_tmp.erase(tsCol_tmp.begin()+i);
					  break ;
					}
				}

			}//================================================================================

		}
		//---------------------------------------------------------------------------------------------
		// // add all tracks that have not been merged :

		debug() << "Still works here" << endmsg;
		for( TrackClusterer::element_vector::iterator it = curSegVec.begin(); it != curSegVec.end() ;++it){

                                // std::cout << "Can I get the track? 1" << std::endl;
			if( (*it)->second == 0 ){

                                // std::cout << "Can I get the track? 2" << std::endl;
                                // std::cout << (*it)->first->edm4hepTrack;
			        edm4hep::Track trk = (*it)->first->edm4hepTrack ;
                                // std::cout << "Can I get the track? 3" << std::endl;
                                // std::cout << trk << std::endl;

				if( copyTrackSegments) {

				        edm4hep::Track t =  trk;

					TrackInfo_of_edm4hepTrack(t) = 0 ; // set extension to 0 to prevent double free ...

					MarTrk_of_edm4hepTrack(t) = 0 ; // dynamic_cast<Track*>( (*it)->first )->ext<MarTrk>() ;

					// FIXME Mingrui debug
					// streamlog_out( DEBUG2 ) << "   create new track from existing LCIO track  - ptr to MarlinTrk : " << t->ext<MarTrk>()  << std::endl ;

					outCol_tmp.push_back( new ClupaPlcioTrack(t) ) ;

				} else {
				        outCol_tmp.push_back( new ClupaPlcioTrack(trk) ) ;
				        debug() << "now outCol_tmp size = " << outCol_tmp.size() << endmsg;
					//remove from segment collection:
					for( int i=0,N=tsCol_tmp.size() ; i<N ; ++i) {
					  if( ( tsCol_tmp.at(i)->edm4hepTrack ) == trk ){
						        // FIXME Do not remove
							// tsCol->removeElementAt( i ) ;
						  //fucd
						  //delete tsCol_tmp[i];
						  tsCol_tmp.erase(tsCol_tmp.begin()+i);
						  break ;
						}
					}
				}


			}
		}
	}

	timer.time( t_merge ) ;

// debug() << "Finish merging works here" << endmsg;

	//---------------------------------------------------------------------------------------------------------
	//    pick up hits from Si trackers
	//---------------------------------------------------------------------------------------------------------

	/*
	   if( _pickUpSiHits ){
	   pickUpSiTrackerHits(outCol) ;
	   }
	   */
	//---------------------------------------------------------------------------------------------------------

	timer.time( t_pickup ) ;

        totalCandidates = outCol_tmp.size();
        debug() << "Total " << totalCandidates << " clupatra tracks to save" << endmsg;
        for (auto t : outCol_tmp) {
	  debug() << "Final track: omega = " << getOmega(t->edm4hepTrack) << endmsg;
	  // omega = getOmega(t->edm4hepTrack);
	  // tree->Fill();
	  outCol->push_back(t->edm4hepTrack);
        }
	for (auto t : tsCol_tmp) {
	  debug() << "Final track segment: omega = " << getOmega(t->edm4hepTrack) << endmsg;
	  tsCol->push_back(t->edm4hepTrack);
	}
	for (auto t : outCol_tmp) {
            delete t;
        }
        for (auto t : tsCol_tmp) {
            delete t;
        }

	_nEvt++ ;

	TrackInfo_of_edm4hepTrack.clear();
	MarTrk_of_edm4hepTrack.clear();
	CluTrk_of_MarTrack.clear();
	MarTrkof.clear();
	GHitof.clear();

	return StatusCode::SUCCESS;
}

// ####### 001
/*************************************************************************************************/

void ClupatraAlg::computeTrackInfo(  edm4hep::ConstTrack lTrk  ){

	if( ! TrackInfo_of_edm4hepTrack(lTrk) )
		TrackInfo_of_edm4hepTrack(lTrk) = new TrackInfoStruct ;

	float r_inner = _gearTPC->getPlaneExtent()[0] ;
	float r_outer = _gearTPC->getPlaneExtent()[1] ;
	float driftLength = _gearTPC->getMaxDriftLength() ;

	// compute z-extend of this track segment
	// auto hv = lTrk->getTrackerHits() ;

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


	// protect against bad tracks
	// if(  tsF == 0 ) return ;
	// if(  tsL == 0 ) return ;
	if (! hasTrackStateAt(lTrk, lcio::TrackState::AtFirstHit)) return; // StatusCode::FAILURE;
	if (! hasTrackStateAt(lTrk, lcio::TrackState::AtFirstHit)) return; // StatueCode::FAILURE;

	edm4hep::TrackState tsF = getTrackStateAt(lTrk, lcio::TrackState::AtFirstHit  ) ;
	edm4hep::TrackState tsL = getTrackStateAt(lTrk, lcio::TrackState::AtLastHit  ) ;

	gear::Vector3D fhPos(tsF.referencePoint[0], tsF.referencePoint[1], tsF.referencePoint[2]) ;
	gear::Vector3D lhPos(tsL.referencePoint[0], tsL.referencePoint[1], tsL.referencePoint[2]) ;


	TrackInfoStruct* ti = TrackInfo_of_edm4hepTrack(lTrk);

	ti->startsInner =  std::abs( fhPos.rho() - r_inner )     <  _trackStartsInnerDist ;        // first hit close to inner field cage
	ti->isCentral   =  std::abs( lhPos.rho() - r_outer )     <  _trackEndsOuterCentralDist ;   // last hit close to outer field cage
	ti->isForward   =  driftLength - std::abs( lhPos.z() )   <  _trackEndsOuterForwardDist  ;  // last hitclose to endcap
	ti->isCurler    =  std::abs( tsF.omega )           >  _trackIsCurlerOmega  ;         // curler segment ( r <~ 1m )

	ti->zMin = zMin ;
	ti->zMax = zMax ;
	ti->zAvg = zAvg ;
	
	//debug() << lTrk.id() << " " << tsF.D0 << " " << tsF.phi << " " << tsF.omega << " " << tsF.Z0 << " " << tsF.tanLambda << endmsg;
}


/*************************************************************************************************/
// void ClupatraAlg::check() {


// check that all Clupatra tracks actually have the four canonical track states set:

// FIXME debug
// streamlog_out( DEBUG5 ) <<   " ------------------- check()  called " << std::endl ;

//  for(  LCIterator<Track> it(   evt, _segmentsOutColName  ) ; EVENT::Track* trk = it.next()  ; ) {
/* FIXME debug
   for(  LCIterator<Track> it(   evt, _outColName  ) ; EVENT::Track* trk = it.next()  ; ) {

   const EVENT::TrackState* ts0 = trk->getTrackState( edm4hep::TrackState::AtIP ) ;
   const EVENT::TrackState* ts1 = trk->getTrackState( edm4hep::TrackState::AtFirstHit ) ;
   const EVENT::TrackState* ts2 = trk->getTrackState( edm4hep::TrackState::AtLastHit ) ;
   const EVENT::TrackState* ts3 = trk->getTrackState( edm4hep::TrackState::AtCalorimeter ) ;


//    streamlog_out( DEBUG2 ) <<  lcshort( trk ) <<  ", " << ts0 <<  ", " << ts1 <<  ", " << ts2 <<  ", " << ts3  << std::endl ;
// FIXME debug
// streamlog_out( DEBUG3 ) <<  " -- " << ts0 <<  ", " << ts1 <<  ", " << ts2 <<  ", " << ts3  << std::endl ;

if( ! ts0 || ! ts1 || ! ts2 || ! ts3 )

// FIXME debug
// streamlog_out( ERROR ) <<  " clupatra track w/ missing track state : " <<  lcshort( trk )
// <<  "  ts0-ts3 : " << ts0 <<  ", " << ts1 <<  ", " << ts2 <<  ", " << ts3  << std::endl ;


}

streamlog_out( DEBUG5 ) <<   " ------------------- check()  done " << std::endl ;
*/

/*************************************************************************************************/
// }



//====================================================================================================

StatusCode ClupatraAlg::finalize(){
  //tree->SaveAs("Tuple.root");

	/*
	 * FIXME debug
	 streamlog_out( MESSAGE )  << "ClupatraProcessor::end()  " << name()
	 << " processed " << _nEvt << " events in " << _nRun << " runs "
	 << std::endl ;
	 */
	return StatusCode::SUCCESS;

}


//====================================================================================================
