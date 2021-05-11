#ifndef clupatra_new_h
#define clupatra_new_h

#include "RuntimeMap.h"
#include "Tracking/TrackingHelper.h"

#include <cmath>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <sstream>
#include <memory>
#include "assert.h"

#include "NNClusterer.h"

#include "lcio.h"
#include "EVENT/TrackerHit.h"
#include "IMPL/TrackImpl.h"
#include "UTIL/Operators.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/ILDConf.h"

#include "gear/GEAR.h"

#include "edm4hep/TrackState.h"

#include "TrackSystemSvc/IMarlinTrack.h"
#include "TrackSystemSvc/IMarlinTrkSystem.h"

// ----- include for verbosity dependend logging ---------
// #include "marlin/VerbosityLevels.h"


/** Helper structs that should go to LCIo to make extraction of layer (and subdetector etc easier )
*/
namespace lcio{
	struct ILDDecoder : public CellIDDecoder<TrackerHit>{
		ILDDecoder() :  lcio::CellIDDecoder<TrackerHit>( lcio::ILDCellID0::encoder_string ) {}
	} ;
	static ILDDecoder ILD_cellID ;

	struct ILDTrackTypeBit{
		static const int SEGMENT ;
		static const int COMPOSITE ;
	} ;
}

namespace clupatra_new{

	/** Small wrapper extension of the LCIO Hit
	*/
	struct ClupaHit {

		ClupaHit() :layer(-1),
		zIndex(-1),
		phiIndex(-1),
		edm4hepHit(0),
		pos(0.,0.,0.) {}
		int layer ;
		int zIndex ;
		int phiIndex ;
		edm4hep::ConstTrackerHit edm4hepHit ;
		gear::Vector3D pos ;

	};

	//  inline lcio::TrackerHit* lcioHit( const ClupaHit* h) { return h->lcioHit ; }


	//------------------ typedefs for elements and clusters ---------

	typedef nnclu::NNClusterer< ClupaHit > Clusterer ;

	typedef Clusterer::element_type Hit ;
	typedef Clusterer::cluster_type CluTrack ;

	typedef Clusterer::element_vector HitVec ;
	typedef Clusterer::cluster_vector CluTrackVec ;

	typedef std::list<Hit*>        HitList ;
	typedef std::vector< HitList > HitListVector ;

	struct ClupaPlcioTrack {
		edm4hep::Track edm4hepTrack ;
                ClupaPlcioTrack(edm4hep::Track edm4hepTrack) : edm4hepTrack(edm4hepTrack) {}
	};

	// typedef GenericHitVec<ClupaHit>      GHitVec ;
	// typedef GenericClusterVec<ClupaHit>  GClusterVec ;

	//------------------------------------------------------------------------------------------

	// struct GHit : lcrtrel::LCExtension<GHit, Hit > {} ;

	//------------------------------------------------------------------------------------------

	// struct MarTrk : lcrtrel::LCExtension<MarTrk, MarlinTrk::IMarlinTrack> {} ;

	//------------------------------------------------------------------------------------------

	// struct DChi2 : lcrtrel::LCFloatExtension<DChi2> {} ;

	//----------------------------------------------------------------

	/** Simple predicate class for computing an index from N bins of the z-coordinate of LCObjects
	 *  that have a float/double* getPostion() method.
	 */
	class ZIndex{
		public:
			/** C'tor takes zmin and zmax - NB index can be negative and larger than N */
			ZIndex( float zmin , float zmax , int n ) : _zmin( zmin ), _zmax( zmax ) , _N (n) {}

			template <class T>
				inline int operator() (T* hit) {

					return (int) std::floor( ( hit->getPosition()[2] - _zmin ) / ( _zmax - _zmin ) * _N ) ;
				}

			inline int index( double z) {  return  (int) std::floor( ( z - _zmin ) / ( _zmax - _zmin ) * _N ) ;  }

		protected:
			ZIndex() {} ;
			float _zmin ;
			float _zmax ;
			int _N ;
	} ;

	//------------------------------------------------------------------------------------------

	struct ZSort {
		inline bool operator()( const Hit* l, const Hit* r) {
			return ( l->first->pos.z() < r->first->pos.z() );
		}
	};


	//------------------------------------------------------------------------------------------

	struct PtSort {  // sort tracks wtr to pt - largest first
		inline bool operator()( const edm4hep::Track l, const edm4hep::Track r) {
			return ( std::abs( getOmega(l) ) < std::abs( getOmega(r) )  );  // pt ~ 1./omega
		}
	};

	//------------------------------------------------------------------------------------------


	/** Add the elements (Hits) from (First,Last) to the HitListVector - vector needs to be initialized, e.g. with
	 *  hLV.resize( MAX_LAYER_INDEX )
	 */
	template <class First, class Last>
		void addToHitListVector( First first, Last last, HitListVector& hLV ){

			while( first != last ){

				// FIXME debug Mingrui
				//streamlog_out( DEBUG ) << "  add hit << "  <<   *first << " to layer " <<  (*first)->first->layer  << std::endl ;

				hLV[ (*first)->first->layer ].push_back( *first )  ;
				++first ;
			}
		}

	//------------------------------------------------------------------------------------------

	/** Predicate class for 'distance' of NN clustering. */
	class HitDistance{
		public:

			HitDistance(float dCut, float caCut = -1.0 ) : _dCutSquared( dCut*dCut ) , _caCut( caCut ) {}

			/** Merge condition: true if distance  is less than dCut */
			inline bool operator()( Hit* h0, Hit* h1){

				if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;

				if( h0->first->layer == h1->first->layer )
					return false ;

				if(  _caCut > 0.  && std::abs( h0->first->layer - h1->first->layer ) == 1 ){

					gear::Vector3D& p0 =  h0->first->pos   ;
					gear::Vector3D& p1 =  h1->first->pos   ;

					double cosAlpha = p0.dot( p1 ) / p0.r() / p1.r()  ;

					// merge hits that seem to come from stiff track from the IP
					//fixme: make parameter
					if( cosAlpha > _caCut ) return true ;
				}

				return ( h0->first->pos - h1->first->pos).r2()  < _dCutSquared ;
			}
		protected:
			HitDistance() ;
			float _dCutSquared, _caCut  ;
	} ;


	// /** Predicate class for 'distance' of NN clustering. */

	// struct HitDistance{  float _dCutSquared ;
	//   HitDistance(float dCut) : _dCutSquared( dCut*dCut ) {}

	//   inline bool operator()( Hit* h0, Hit* h1){

	//     if( h0->first->layer == h1->first->layer )
	// 	return false ;

	//     return ( h0->first->pos - h1->first->pos).r2()
	// 	< _dCutSquared ;
	//   }
	// } ;



	//-----------------------------------------------

	struct PLCIOTrackConverter{

		bool UsePropagate ;
		PLCIOTrackConverter() : UsePropagate(false ) {}

		edm4hep::Track operator() (CluTrack* c) ;

	} ;

	//------------------------------------------------------------------------------------------
	/** Predicate class for identifying small clusters. */
	struct ClusterSize {
		ClusterSize( unsigned n) : _n(n) {}
		bool operator()(const CluTrack* cl) const {   return cl->size() < _n ;  }
		unsigned _n ;
	};

	//------------------------------------------------------------------------------------------
	/** Predicate class for identifying clusters with duplicate pad rows - returns true
	 *  if the fraction of duplicate hits is larger than 'fraction'.
	 */
	struct DuplicatePadRows{

		unsigned _N ;
		float _f ;
		DuplicatePadRows(unsigned nPadRows, float fraction) : _N( nPadRows), _f( fraction )  {}

		bool operator()(const CluTrack* cl) const {

			// check for duplicate layer numbers
			std::vector<int> hLayer( _N )  ;
			typedef CluTrack::const_iterator IT ;

			unsigned nHit = 0 ;
			for(IT it=cl->begin() ; it != cl->end() ; ++it ) {
				++ hLayer[ (*it)->first->layer]   ;
				++ nHit ;
			}
			unsigned nDuplicate = 0 ;
			for(unsigned i=0 ; i < _N ; ++i ) {
				if( hLayer[i] > 1 )
					nDuplicate += hLayer[i] ;
			}
			return double(nDuplicate)/nHit > _f ;
		}
	};

	//------------------------------------------------------------------------------------------

	// helper for sorting cluster wrt layerID
	struct LayerSortOut{
		bool operator()( const Hit* l, const Hit* r) { return l->first->layer < r->first->layer ; }
	} ;
	struct LayerSortIn{
		bool operator()( const Hit* l, const Hit* r) { return l->first->layer > r->first->layer ; }
	} ;




	//------------------------------------------------------------------------------------------

	struct IMarlinTrkFitter{

		MarlinTrk::IMarlinTrkSystem* _ts ;
		double _maxChi2Increment ;

		IMarlinTrkFitter(MarlinTrk::IMarlinTrkSystem* ts, double maxChi2Increment=DBL_MAX ) :
			_ts( ts ) ,
			_maxChi2Increment(maxChi2Increment) {}


		MarlinTrk::IMarlinTrack* operator() (CluTrack* clu) ;
	};

	//-------------------------------------------------------------------------------------

	/** Try to add hits from hLV (hit lists per layer) to the cluster. The cluster needs to have a fitted KalTrack associated to it.
	 *  Hits are added if the resulting delta Chi2 is less than dChiMax - a maxStep is the maximum number of steps (layers) w/o
	 *  successfully merging a hit.
	 */
	int addHitsAndFilter( CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut, unsigned maxStep, ZIndex& zIndex,  bool backward=false,
			MarlinTrk::IMarlinTrkSystem* trkSys=0) ;
	//------------------------------------------------------------------------------------------

	/** Try to add a hit from the given HitList in layer of subdetector to the track.
	 *  A hit is added if the resulting delta Chi2 is less than dChiMax.
	 */
	bool addHitAndFilter( int detectorID, int layer, CluTrack* clu, HitListVector& hLV , double dChiMax, double chi2Cut) ;

	//------------------------------------------------------------------------------------------
	/** Split up clusters that have a hit multiplicity of 2,3,4,...,N in at least layersWithMultiplicity.
	*/
	void split_multiplicity( Clusterer::cluster_list& cluList, int layersWithMultiplicity , int N=5) ;

	//------------------------------------------------------------------------------------------
	/** Returns the number of rows where cluster clu has i hits in mult[i] for i=1,2,3,4,.... -
	 *  mult[0] counts all rows that have hits
	 */
	void getHitMultiplicities( CluTrack* clu, std::vector<int>& mult ) ;

	//------------------------------------------------------------------------------------------

	/** Split the cluster into two clusters.
	*/
	void create_two_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ) ;

	//------------------------------------------------------------------------------------------

	/** Split the cluster into three clusters.
	*/
	void create_three_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ) ;


	/** Split the cluster into N clusters.
	*/
	void create_n_clusters( Clusterer::cluster_type& clu, Clusterer::cluster_list& cluVec ,  unsigned n ) ;



	//=======================================================================================

	struct TrackInfoStruct{
		TrackInfoStruct() : zMin(0.), zAvg(0.), zMax(0.), startsInner(false), isCentral(false), isForward(false), isCurler(false) {}
		float zMin ;
		float zAvg ;
		float zMax ;
		bool startsInner ;
		bool isCentral   ;
		bool isForward   ;
		bool isCurler    ;
	} ;
	// struct TrackInfo : lcrtrel::LCOwnedExtension<TrackInfo, TrackInfoStruct> {} ;
	typedef TrackInfoStruct TrackInfo;
	/** Helper class to compute track segment properties.
	*/
	struct ComputeTrackerInfo{
		void operator()( edm4hep::Track o );
	};

	//=======================================================================================

	/** helper class for merging split track segments */

	class TrackSegmentMerger{

		public:
			/** C'tor takes merge distance */
			TrackSegmentMerger(float chi2Max,  MarlinTrk::IMarlinTrkSystem* trksystem, float b) : _chi2Max( chi2Max ) , _trksystem( trksystem), _b(b) {}

			float _chi2Max ;
			MarlinTrk::IMarlinTrkSystem* _trksystem ;
			float _b ;

			/** Merge condition: ... */
			inline bool operator()( nnclu::Element<ClupaPlcioTrack>* h0, nnclu::Element<ClupaPlcioTrack>* h1){

				edm4hep::Track trk0 = h0->first->edm4hepTrack ;
				edm4hep::Track trk1 = h1->first->edm4hepTrack ;


				// protect against merging multiple segments (and thus complete tracks)
				if(  h0->second || h1->second )
					return false ;

				// const TrackInfoStruct* ti0 =  trk0->ext<TrackInfo>() ;
				// const TrackInfoStruct* ti1 =  trk1->ext<TrackInfo>() ;

				// const lcio::TrackState* tsF0 = trk0->getTrackState( lcio::TrackState::AtFirstHit  ) ;
				// const lcio::TrackState* tsL0 = trk0->getTrackState( lcio::TrackState::AtLastHit  ) ;

				// const lcio::TrackState* tsF1 = trk1->getTrackState( lcio::TrackState::AtFirstHit  ) ;
				// const lcio::TrackState* tsL1 = trk1->getTrackState( lcio::TrackState::AtLastHit  ) ;

				// --- get three hits per track : first last midddle

				unsigned nhit0 = trk0.trackerHits_size() ;
				unsigned nhit1 = trk1.trackerHits_size() ;

				edm4hep::ConstTrackerHit thf0 = trk0.getTrackerHits( 0 ) ;
				edm4hep::ConstTrackerHit thf1 = trk1.getTrackerHits( 0 ) ;

				edm4hep::ConstTrackerHit thl0 = trk0.getTrackerHits( nhit0 - 1 ) ;
				edm4hep::ConstTrackerHit thl1 = trk1.getTrackerHits( nhit1 - 1 ) ;

				// lcio::TrackerHit* thm1 = trk1->getTrackerHits()[ nhit1 / 2 ] ;
				// lcio::TrackerHit* thm0 = trk0->getTrackerHits()[ nhit0 / 2 ] ;

				// int lthf0 = ILD_cellID(  thf0 )[ lcio::ILDCellID0::layer ] ;
				// int lthf1 = ILD_cellID(  thf1 )[ lcio::ILDCellID0::layer ] ;

				// int lthl0 = ILD_cellID(  thl0 )[ lcio::ILDCellID0::layer ] ;
				// int lthl1 = ILD_cellID(  thl1 )[ lcio::ILDCellID0::layer ] ;
				int lthf0 = getLayer(thf0);
				int lthf1 = getLayer(thf1);
				int lthl0 = getLayer(thl0);
				int lthl1 = getLayer(thl1);

				//      if( lthf0 <= lthl1 && lthf1 <= lthl0 )   return false ;

				// allow the track segements to overlap slightly  - FIXME: make a parameter ...
				const int overlapRows = 4 ;
				if( lthf0 + overlapRows <= lthl1 && lthf1  + overlapRows  <= lthl0 )   return false ;

				// now we take the larger segment and see if we can add the three hits from the other segment...

				edm4hep::Track trk = ( nhit0 > nhit1 ? trk0 :  trk1 ) ;
				edm4hep::Track oth = ( nhit0 > nhit1 ? trk1 :  trk0 ) ;

				bool  outward = ( nhit0 > nhit1  ?  lthl0 <= lthf1 + overlapRows :  lthl1 <= lthf0 + overlapRows ) ;

				unsigned n = oth.trackerHits_size() ;

				edm4hep::ConstTrackerHit th0 =  ( outward ? oth.getTrackerHits(0) : oth.getTrackerHits(n-1) ) ;
				edm4hep::ConstTrackerHit th1 =              (oth.getTrackerHits(n/2) );
				edm4hep::ConstTrackerHit th2 =  ( outward ? oth.getTrackerHits(n-1) : oth.getTrackerHits(0) );


				// track state at last hit migyt be rubish....
				//      const lcio::TrackState* ts = ( outward ? trk->getTrackState( lcio::TrackState::AtLastHit  ) : trk->getTrackState( lcio::TrackState::AtFirstHit  ) ) ;
				const edm4hep::TrackState ts = getTrackStateAt(trk, lcio::TrackState::AtFirstHit  )  ;



				// FIXME Mingrui debug
				// streamlog_out( DEBUG3 ) << " *******  TrackSegmentMerger : will extrapolate track " << ( outward ? " outwards\t" : " inwards\t" )
				//<<  lcio::lcshort( trk  ) << "     vs:  [" <<   std::hex << oth->id() << std::dec << "]"  << std::endl ;

				// if( trk->id() == 0x0004534d &&  oth->id() == 0x000454a6 ){
				// 	streamlog_out( DEBUG3 )  << " &&&&&&&&&&&&&& Track 1 : \n" << *trk
				// 				 << " &&&&&&&&&&&&&& Track 2 : \n" << *oth
				// 				 <<  std::endl ;
				// }

				std::auto_ptr<MarlinTrk::IMarlinTrack> mTrk( _trksystem->createTrack()  ) ;

				int nHit = trk.trackerHits_size() ;

				if( nHit == 0 /*|| ts ==0*/ ) // FIXME Mingrui Since it is no pointer
					return false ;

				// float initial_chi2 = trk->getChi2() ;
				// float initial_ndf  = trk->getNdf() ;

				// FIXME Mingrui debug
				// streamlog_out( DEBUG3  )  << "               -- extrapolate TrackState : " << lcshort( ts )    << std::endl ;

				edm4hep::ConstTrackerHit ht = trk.getTrackerHits(0); 
				//need to add a dummy hit to the track
				mTrk->addHit( ht ) ;  // is this the right hit ??????????

				mTrk->initialise( ts ,  _b ,  MarlinTrk::IMarlinTrack::backward ) ;

				double deltaChi ;
				int addHit = 0 ;

				//-----   now try to add the three hits : ----------------
				addHit = mTrk->addAndFit(  th0 , deltaChi, _chi2Max ) ;

				/*
				 * FIXME Mingrui debug
				 streamlog_out( DEBUG3 ) << "    ****  adding first hit : " <<  gear::Vector3D( th0->getPosition() )
				 << "         added : " << MarlinTrk::errorCode( addHit )
				 << "         deltaChi2: " << deltaChi
				 << std::endl ;
				 */

				if( addHit !=  MarlinTrk::IMarlinTrack::success ) return false ;

				//---------------------
				addHit = mTrk->addAndFit(  th1 , deltaChi, _chi2Max ) ;

				/* FIXME Mingrui debug
				   streamlog_out( DEBUG3 ) << "    ****  adding second hit : " <<  gear::Vector3D( th1->getPosition() )
				   << "         added : " << MarlinTrk::errorCode( addHit )
				   << "         deltaChi2: " << deltaChi
				   << std::endl ;
				   */

				if( addHit !=  MarlinTrk::IMarlinTrack::success ) return false ;

				//--------------------
				addHit = mTrk->addAndFit(  th2 , deltaChi, _chi2Max ) ;

				/* FIXME Mingrui debug
				   streamlog_out( DEBUG3 ) << "    ****  adding third hit : " <<  gear::Vector3D( th2->getPosition() )
				   << "         added : " << MarlinTrk::errorCode( addHit )
				   << "         deltaChi2: " << deltaChi
				   << std::endl ;
				   */

				if( addHit !=  MarlinTrk::IMarlinTrack::success ) return false ;



				////////////
				return true ;
				///////////
			}

	};
	//=======================================================================================

	/** helper class for merging track segments, based on circle (and tan lambda) */

	class TrackCircleDistance{

		public:
			/** C'tor takes merge distance */
			TrackCircleDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut){}

			/** Merge condition: ... */
			bool operator() ( nnclu::Element<ClupaPlcioTrack>* h0, nnclu::Element<ClupaPlcioTrack>* h1);
  protected:
    float _dCutSquared ;
    float _dCut ;
			//=======================================================================================
	};

	struct TrackZSort {  // sort tracks wtr to abs(z_average )
		bool operator()( edm4hep::Track l, edm4hep::Track r);
	};


	//=======================================================================================

	/** Helper class that creates an Elements for an LCOjects of type T.
	*/
	template <class T>
		struct MakePLCIOElement{
			nnclu::Element<T>*  operator()( T* o) { return new nnclu::Element<T>(o) ;  }
		} ;


	//=======================================================================================

	class Timer{
		public:
			Timer(){
				_clocks.reserve( 100 ) ;
				_names.reserve( 100 ) ;

				_clocks.push_back(0) ;
				_names.push_back(  "start"  ) ;
			}
			unsigned registerTimer( const std::string& name ){
				_clocks.push_back(0) ;
				_names.push_back( name ) ;
				return _clocks.size() - 1 ;
			}

			void time(unsigned index){
				_clocks[ index ] = clock() ;
			}
			void start() { time(0) ; }


			std::string toString(){

				std::stringstream s ;

				s << " ============= Timer ================================ "  << std::endl ;
				unsigned N=_clocks.size() ;
				for( unsigned i=1 ;  i < N ; ++i){
					s << "    " << _names[i] << " : " <<  double(_clocks[i] - _clocks[i-1] ) / double(CLOCKS_PER_SEC) << std::endl ;
				}
				s << "         Total  : " <<  double(_clocks[N-1] - _clocks[0] ) / double(CLOCKS_PER_SEC) << std::endl ;
				s << " ==================================================== "  << std::endl ;

				return s.str() ;
			}
		protected:
			std::vector< clock_t> _clocks ;
			std::vector< std::string > _names ;
	};


}
#endif
