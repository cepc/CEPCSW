#ifndef VXDTrack_h
#define VXDTrack_h

#include "edm4hep/Track.h"
#include "edm4hep/TrackerHit.h"
#include "TrackSystemSvc/IMarlinTrkSystem.h"
#include "TrackSystemSvc/IMarlinTrack.h"

#include <vector>

#include "ILDImpl/IVXDHit.h"
#include "ILDImpl/IMiniVector.h"
#include "KiTrack/ITrack.h"

#include "Tools/Fitter.h"

//#include "SpacePointBuilder.h"
// CLHEP tools
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"


namespace KiTrackMarlin{
  /** A class for ITracks containing an lcio::Track at core
   */
  class VXDTrack : public ITrack {
  public:
      
    /** @param trkSystem An IMarlinTrkSystem, which is needed for fitting of the tracks
     */
    VXDTrack( MarlinTrk::IMarlinTrkSystem* trkSystem );
    
    /** @param hits The hits the track consists of 
     * @param trkSystem An IMarlinTrkSystem, which is needed for fitting of the tracks
     */
    //VXDTrack( std::vector< IVXDHit* > hits , MarlinTrk::IMarlinTrkSystem* trkSystem );
    VXDTrack( std::vector< IMiniVector* > hits , MarlinTrk::IMarlinTrkSystem* trkSystem );
    VXDTrack( const VXDTrack& f );
    VXDTrack & operator= (const VXDTrack & f);
          
    /** @return a track in the lcio format
     */
    edm4hep::Track* getLcioTrack(){ return ( _lcioTrack );}
    
    //void addHit( IVXDHit* hit );
    void addHit( IMiniVector* MV );
    
    virtual double getNdf() const { return _lcioTrack->getNdf(); }
    virtual double getChi2() const { return _lcioTrack->getChi2(); }
    virtual double getChi2Prob() const { return _chi2Prob; }

    virtual std::vector< IHit* > getHits() const {
      std::vector<IHit*> hits; 
      for(unsigned i=0; i<_hits.size();i++) hits.push_back( _hits[i] ); 
      return hits;
    }
 
    virtual std::vector< IMiniVector* > getMVs() const {
      std::vector<IMiniVector*> mvhits; 
      for(unsigned i=0; i<_hits.size();i++) mvhits.push_back( _hits[i] ); 
      return mvhits;
    }
     
    virtual double getQI() const;
    
    /** Fits the track and sets chi2, Ndf etc.
     */
    virtual void fit() ;
    
    virtual ~VXDTrack(){ delete _lcioTrack; }
      
  protected:
    /** the hits the track consists of
     */
    //std::vector< IVXDHit* > _hits;
    std::vector< IMiniVector* > _hits;     
      
    edm4hep::Track* _lcioTrack;
      
    // for fitting
    MarlinTrk::IMarlinTrkSystem* _trkSystem;
        
    double _chi2Prob;
  };
}
#endif


