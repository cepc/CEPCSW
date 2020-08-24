#include "TrackSystemSvc/IMarlinTrkSystem.h"
#include <sstream>

namespace MarlinTrk{
  
  void IMarlinTrkSystem::setOption(unsigned CFGOption, bool val) {
    _cfg.setOption( CFGOption, val ) ;
  }
  
  
  bool IMarlinTrkSystem::getOption( unsigned CFGOption) {
    return _cfg[ CFGOption] ;
  }
  
  
  std::string IMarlinTrkSystem::getOptions() {
    
    std::stringstream ss ;
    ss << _cfg ;
    return ss.str() ;
  }
  
  void IMarlinTrkSystem::registerOptions() {
    
    _cfg.registerOption( IMarlinTrkSystem::CFG::useQMS,  "useMultipleScattering", true) ;
    _cfg.registerOption( IMarlinTrkSystem::CFG::usedEdx, "useEnergyLoss", true) ;
    _cfg.registerOption( IMarlinTrkSystem::CFG::useSmoothing, "useSmoothingInFit", false) ;
    
    
  }
  
  
}
