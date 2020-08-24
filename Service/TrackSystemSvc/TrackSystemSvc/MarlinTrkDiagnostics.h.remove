#ifndef MarlinTrkDiagnostics_h
#define MarlinTrkDiagnostics_h 

//// switch to turn on diagnostics code
//#define MARLINTRK_DIAGNOSTICS_ON 

#ifdef MARLINTRK_DIAGNOSTICS_ON

#include "lcio.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TrackerHit.h"


namespace MarlinTrk{

  
  // LCIO Extension creating a pointer to the simhit for trackerhits 
  struct MCTruth4HitExtStruct{
    MCTruth4HitExtStruct() : simhit(0) {}
    EVENT::SimTrackerHit* simhit;
  } ; 
  struct MCTruth4HitExt : lcio::LCOwnedExtension<MCTruth4HitExt, MCTruth4HitExtStruct> {} ;
  
  // fills a vector of MCParticle pointers with the MCParticles assosiated with the provided tracker hit using MCTruth4HitExtStruct
  void getMCParticlesForTrackerHit(EVENT::TrackerHit* trkhit, std::vector<EVENT::MCParticle*>& mcps) ;
  
}

#endif

#endif
