#ifndef DiagnosticsController_h
#define DiagnosticsController_h 

#include "MarlinTrk/MarlinTrkDiagnostics.h"

#ifdef MARLINTRK_DIAGNOSTICS_ON

class TFile;
class TH1F;
class TTree;
class TKalMatrix;

namespace EVENT {
  class MCParticle;
}

class ILDVTrackHit;
class TKalTrackSite;
class MarlinTrkNtuple;

namespace MarlinTrk{

  class MarlinKalTestTrack;
  class IMarlinTrack;
 
  
  class DiagnosticsController {
    
  public:
    
    /** constructor */
    DiagnosticsController(); 
    
    /** Destructor */
    virtual ~DiagnosticsController();  
    
    
    void init(std::string root_file_name, std::string root_Tree_name, bool _recording_on=true ) ;
    
    void new_track(MarlinKalTestTrack* trk) ;
    
    void set_intial_track_parameters(double d0, double phi0, double omega, double z0, double tanL, double pivot_x, double pivot_y, double pivot_z, TKalMatrix& cov);
    
    void record_site(ILDVTrackHit* hit, TKalTrackSite* site);
    
    void record_rejected_site(ILDVTrackHit* hit, TKalTrackSite* site);
    
    void skip_current_track(); 
      
    void end_track() ;
    
    void end();
    
        
  private:
           
    DiagnosticsController(const DiagnosticsController&) ;                 // Prevent copy-construction
    DiagnosticsController& operator=(const DiagnosticsController&) ;      // Prevent assignment
    
    void clear_track_record();
    
    bool _initialised;
    bool _recording_on;
    
    int _ntracks_written;
    int _ntracks_skipped;
    
    std::string _root_file_name;
    std::string _root_tree_name;

    TFile* _root_file;
    TTree* _tree;
    MarlinTrkNtuple* _track_record;
    
    MarlinKalTestTrack* _current_track;
    
    EVENT::MCParticle* _currentMCP;

    bool _mcpInfoStored;
    bool _skip_track;
    
  
  };


}

#endif

#endif
