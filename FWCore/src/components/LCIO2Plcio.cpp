#include "LCIO2Plcio.h"

#include "lcio.h"
#include "EVENT/LCCollection.h" 
#include "plcio/MCParticle.h"
#include "plcio/MCParticleConst.h"
#include "plcio/MCParticleCollection.h"
#include "plcio/LCRunHeader.h"
#include "plcio/LCRunHeaderCollection.h"
#include "EVENT/MCParticle.h"

typedef std::vector<EVENT::MCParticle*> MCParticleVec ;

//std::map<std::string, fptr> LCIO2Plcio::map_cvt;
CollectionsMap LCIO2Plcio::map_cols;

void lp_info(std::string s){
  printf("[LCIO2Plcio]:\t\t%s.\n", &s[0]);
}

LCIO2Plcio::LCIO2Plcio(){
  map_cvt.insert(std::make_pair<std::string, fptr>("MCParticle", Convertor_MCParticle));
  map_cvt.insert(std::make_pair<std::string, fptr>("LCRunHeader", Convertor_LCRunHeader));
  map_cvt.insert(std::make_pair<std::string, fptr>("SimTrackerHit", Convertor_SimTrackerHit));
  map_cvt.insert(std::make_pair<std::string, fptr>("SimCalorimeterHit", Convertor_SimCalorimeterHit));
}
LCIO2Plcio::LCIO2Plcio(EVENT::LCCollection* collection){
  LCIO2Plcio();
}

podio::CollectionBase* LCIO2Plcio::Convertor_LCRunHeader(EVENT::LCCollection* lc_col){
  plcio::LCRunHeaderCollection* pl_col = new plcio::LCRunHeaderCollection();
//  plcio::LCRunHeader pl_rh = (plcio::LCRunHeader) pl_col->create();
  return pl_col;
}

void LCIO2Plcio::void_Core_MCParticle(EVENT::MCParticle* lc_mcp, plcio::MCParticle& pl_mcp){

  pl_mcp.setPDG( lc_mcp->getPDG() );
  pl_mcp.setGeneratorStatus( lc_mcp->getGeneratorStatus() );
  pl_mcp.setSimulatorStatus( lc_mcp->getSimulatorStatus() );
  pl_mcp.setCharge( lc_mcp->getCharge() );
  pl_mcp.setTime( lc_mcp->getTime() );
  pl_mcp.setMass( lc_mcp->getMass() );
  pl_mcp.setStopped( lc_mcp->isStopped() );
  pl_mcp.setOverlay( lc_mcp->isOverlay() );
  pl_mcp.setBackscatter( lc_mcp->isBackscatter() );
  pl_mcp.setDecayedInTracker( lc_mcp->isDecayedInTracker() );
  pl_mcp.setDecayedInCalorimeter( lc_mcp->isDecayedInCalorimeter() );
  pl_mcp.setCreatedInSimulation( lc_mcp->isCreatedInSimulation() );
  pl_mcp.setVertexIsNotEndpointOfParent( lc_mcp->vertexIsNotEndpointOfParent() );
  pl_mcp.setHasLeftDetector( lc_mcp->hasLeftDetector() );

  pl_mcp.setSpin( plcio::FloatThree( lc_mcp->getSpin() ) );
  pl_mcp.setColorFlow( plcio::IntTwo( lc_mcp->getColorFlow() ) );
  pl_mcp.setVertex( plcio::DoubleThree( lc_mcp->getVertex()));
  pl_mcp.setEndpoint( plcio::DoubleThree( lc_mcp->getEndpoint() ) );

  // send float ptr as parameter to FloatThree.
  float plcio_m[3];
  const double* lcio_m = lc_mcp->getMomentum();
  plcio_m[0] = lcio_m[0];
  plcio_m[1] = lcio_m[1];
  plcio_m[2] = lcio_m[2];
  pl_mcp.setMomentum( plcio::FloatThree( plcio_m ) );

  lcio_m = lc_mcp->getMomentumAtEndpoint();
  plcio_m[0] = lcio_m[0];
  plcio_m[1] = lcio_m[1];
  plcio_m[2] = lcio_m[2];
  pl_mcp.setMomentumAtEndpoint( plcio::FloatThree( plcio_m ) );
}

podio::CollectionBase* LCIO2Plcio::Convertor_MCParticle(EVENT::LCCollection* lc_col){
  return Core_MCParticle(lc_col);
}

plcio::MCParticleCollection* LCIO2Plcio::Core_MCParticle(EVENT::LCCollection* lc_col){
  plcio::MCParticleCollection* pl_col = new plcio::MCParticleCollection();

  // Convert basic info from LCIO to plcio;
  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
      EVENT::MCParticle* lc_mcp = (EVENT::MCParticle*) lc_col->getElementAt(i) ;
      plcio::MCParticle pl_mcp = (plcio::MCParticle) pl_col->create();

      void_Core_MCParticle(lc_mcp, pl_mcp);
      // dealing each of 'Parents' of lcio::MCParticle ;
      const MCParticleVec& veclc = lc_mcp->getParents();
      for(unsigned j=0; j<veclc.size(); j++){
        EVENT::MCParticle* vlcreg = veclc[j];
        for(unsigned k=0; k<i; ++k){
          if(((EVENT::MCParticle*) lc_col->getElementAt(k)) == vlcreg){
            // A loop for plcio's MCParticleCollection to recover plcio's relationship;
            plcio::MCParticle plprt = pl_col->at(k);
            pl_mcp.addParent(plprt);
            plprt.addDaughter(pl_mcp);        
	  }
	}
      }
  }
  return pl_col;
}

//void LCIO2Plcio::setPlcioMCParticleCollection(plcio::MCParticleCollection* pl_col){
//  hitcol_pl = pl_col;
//}

//void LCIO2Plcio::setLCIOMCParticleCollection(EVENT::LCCollection* lc_col){
//  hitcol_lc = lc_col;
//}


podio::CollectionBase* LCIO2Plcio::Convertor_SimTrackerHit(EVENT::LCCollection* lc_col){
  plcio::SimTrackerHitCollection* pl_col = new plcio::SimTrackerHitCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::SimTrackerHit* lc_sth = (EVENT::SimTrackerHit*) lc_col->getElementAt(i) ;
    plcio::SimTrackerHit pl_sth = (plcio::SimTrackerHit) pl_col->create();

    pl_sth.setCellID0( lc_sth->getCellID0() );
    pl_sth.setCellID1( lc_sth->getCellID1() );
    pl_sth.setEDep( lc_sth->getEDep() );
    pl_sth.setTime( lc_sth->getTime() );
    pl_sth.setPathLength( lc_sth->getPathLength() );
    pl_sth.setQuality( lc_sth->getQuality() );
    pl_sth.setPosition( lc_sth->getPosition() );
    pl_sth.setMomentum( lc_sth->getMomentum() );
    pl_sth.setOverlay( lc_sth->isOverlay() );
    pl_sth.setProducedBySecondary( lc_sth->isProducedBySecondary() );

    // Looping the LCIO::MCParticleCollection to pick Particle for Hits;
    CollectionsVec vec_mcp;
    vec_mcp = map_cols["MCParticle"];
    EVENT::LCCollection* hitcol_lc = (EVENT::LCCollection*) vec_mcp[0].first;
    plcio::MCParticleCollection* hitcol_pl = (plcio::MCParticleCollection*) vec_mcp[0].second;

    int index = -1;
    // search corresponding MCParticleCollection*;
    EVENT::MCParticle* mcptr_reg = lc_sth->getMCParticle();
    for( unsigned j=0, M=hitcol_lc->getNumberOfElements(); j<M; ++j){
      EVENT::MCParticle* mcpin_lc = (EVENT::MCParticle*) hitcol_lc->getElementAt(j);
      if( mcpin_lc == mcptr_reg ){
	index = j;
      }
    }
//    if( index != -1) lp_info("Cedar: Convertor SimTrackerHit Success.");
    bool is_empty = false;
    if( index == -1){
      if(mcptr_reg == nullptr) is_empty = true;
//      lp_info("Convertor SimTrackerHit Problem.");
//      printf("LCCollection size: %d\n", hitcol_lc->getNumberOfElements());
//      printf("MCParticleCollection size: %d\n", hitcol_pl->size());
//      printf("Index: %d\n", index);
    }
    if( is_empty == false )
      pl_sth.setMCParticle( hitcol_pl->at(index) );

  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_SimCalorimeterHit(EVENT::LCCollection* lc_col){
  plcio::SimCalorimeterHitCollection* pl_col = new plcio::SimCalorimeterHitCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::SimCalorimeterHit* lc_sch = (EVENT::SimCalorimeterHit*) lc_col->getElementAt(i);
    plcio::SimCalorimeterHit pl_sch = (plcio::SimCalorimeterHit) pl_col->create();

    pl_sch.setCellID0( lc_sch->getCellID0() );
    pl_sch.setCellID1( lc_sch->getCellID1() );
    pl_sch.setEnergy( lc_sch->getEnergy() );

    float plcio_p[3];
    const float* lcio_p = lc_sch->getPosition();
    plcio_p[0] = lcio_p[0];
    plcio_p[1] = lcio_p[1];
    plcio_p[2] = lcio_p[2];
    pl_sch.setPosition( plcio::FloatThree( plcio_p ) );

    // converting from lc_sch to pl_sch on the contribution variables;
    for( unsigned j=0, N=lc_sch->getNMCContributions(); j<N; j++){
      plcio_p[0] = lc_sch->getStepPosition(j)[0];
      plcio_p[1] = lc_sch->getStepPosition(j)[1];
      plcio_p[2] = lc_sch->getStepPosition(j)[2];

      plcio::ConstCaloHitContribution tmp(
	lc_sch->getPDGCont(j), lc_sch->getEnergyCont(j),
	lc_sch->getTimeCont(j), plcio::FloatThree(plcio_p) 
      );
      pl_sch.addContribution( tmp );
    }
  }

  return pl_col;
}

bool LCIO2Plcio::isReady(const std::string& TypeName){
  if( TypeName == "SimTrackerHit" ){
    std::vector<std::string>::iterator it = find(
                  vec_Types.begin(), vec_Types.end(), "MCParticle"
    );
    if( it != vec_Types.end() )
      return true;
  }
  return false;
}

podio::CollectionBase* LCIO2Plcio::Convertor_getPlcio(EVENT::LCCollection* lc_col){

  podio::CollectionBase* collection(nullptr);
  TypeName = lc_col->getTypeName();

//  lp_info("Converting "+TypeName);
  bool ready = true; 
  if( !ready ){
    lp_info("Not ready yet.");
    lp_info("Please put MCParticle firstly in the python file, Please");
    return nullptr;
  }

  fptr fp;
  if( map_cvt.find(TypeName) == map_cvt.end()){
    lp_info("unrecognized "+TypeName);
  }else{
    fp = *map_cvt[TypeName];
  }

  fp = *map_cvt[TypeName];
  collection = fp(lc_col);

// maintain map<TypeName, CollVec>;
    map_cols[TypeName].push_back(
      std::pair<EVENT::LCCollection*, podio::CollectionBase*>(lc_col, collection)
    );

//  lp_info("done.");
  return collection;
}
