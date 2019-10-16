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

  // Convert basic info from LCIO to plcio;
  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
      EVENT::LCRunHeader* lc_var = (EVENT::LCRunHeader*) lc_col->getElementAt(i) ;
      plcio::LCRunHeader pl_var = (plcio::LCRunHeader) pl_col->create();
      pl_var.setRunNumber( lc_var->getRunNumber() );
      pl_var.setDetectorName( lc_var->getDetectorName() );
      pl_var.setDescription( lc_var->getDescription() );
      
      std::vector<std::string> vec_dct = *(lc_var->getActiveSubdetectors());
      for( unsigned j=0,N=vec_dct.size(); j<N; j++){
        pl_var.addActiveSubdetector( vec_dct[j] );
      }
  }
  return pl_col;
}

void LCIO2Plcio::setMCParticle(EVENT::MCParticle* lc_var, plcio::MCParticle& pl_var){

  pl_var.setPDG( lc_var->getPDG() );
  pl_var.setGeneratorStatus( lc_var->getGeneratorStatus() );
  pl_var.setSimulatorStatus( lc_var->getSimulatorStatus() );
  pl_var.setCharge( lc_var->getCharge() );
  pl_var.setTime( lc_var->getTime() );
  pl_var.setMass( lc_var->getMass() );
  pl_var.setStopped( lc_var->isStopped() );
  pl_var.setOverlay( lc_var->isOverlay() );
  pl_var.setBackscatter( lc_var->isBackscatter() );
  pl_var.setDecayedInTracker( lc_var->isDecayedInTracker() );
  pl_var.setDecayedInCalorimeter( lc_var->isDecayedInCalorimeter() );
  pl_var.setCreatedInSimulation( lc_var->isCreatedInSimulation() );
  pl_var.setVertexIsNotEndpointOfParent( lc_var->vertexIsNotEndpointOfParent() );
  pl_var.setHasLeftDetector( lc_var->hasLeftDetector() );

  pl_var.setSpin( plcio::FloatThree( lc_var->getSpin() ) );
  pl_var.setColorFlow( plcio::IntTwo( lc_var->getColorFlow() ) );
  pl_var.setVertex( plcio::DoubleThree( lc_var->getVertex()));
  pl_var.setEndpoint( plcio::DoubleThree( lc_var->getEndpoint() ) );

  // send float ptr as parameter to FloatThree.
  float plcio_v[3];
  const double* lcio_v = lc_var->getMomentum();
  plcio_v[0] = lcio_v[0];
  plcio_v[1] = lcio_v[1];
  plcio_v[2] = lcio_v[2];
  pl_var.setMomentum( plcio::FloatThree( plcio_v ) );

  lcio_v = lc_var->getMomentumAtEndpoint();
  plcio_v[0] = lcio_v[0];
  plcio_v[1] = lcio_v[1];
  plcio_v[2] = lcio_v[2];
  pl_var.setMomentumAtEndpoint( plcio::FloatThree( plcio_v ) );
}

podio::CollectionBase* LCIO2Plcio::Convertor_MCParticle(EVENT::LCCollection* lc_col){
  plcio::MCParticleCollection* pl_col = new plcio::MCParticleCollection();

  // Convert basic info from LCIO to plcio;
  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
      EVENT::MCParticle* lc_var = (EVENT::MCParticle*) lc_col->getElementAt(i) ;
      plcio::MCParticle pl_var = (plcio::MCParticle) pl_col->create();

      setMCParticle(lc_var, pl_var);
      // dealing each of 'Parents' of lcio::MCParticle ;
      const MCParticleVec& veclc = lc_var->getParents();
      for(unsigned j=0; j<veclc.size(); j++){
        EVENT::MCParticle* vlcreg = veclc[j];
        for(unsigned k=0; k<i; ++k){
          if(((EVENT::MCParticle*) lc_col->getElementAt(k)) == vlcreg){
            // A loop for plcio's MCParticleCollection to recover plcio's relationship;
            plcio::MCParticle mcprt = pl_col->at(k);
            pl_var.addParent(mcprt);
            mcprt.addDaughter(pl_var);        
	  }
	}
      }
  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_SimTrackerHit(EVENT::LCCollection* lc_col){
  plcio::SimTrackerHitCollection* pl_col = new plcio::SimTrackerHitCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::SimTrackerHit* lc_var = (EVENT::SimTrackerHit*) lc_col->getElementAt(i) ;
    plcio::SimTrackerHit pl_var = (plcio::SimTrackerHit) pl_col->create();

    pl_var.setCellID0( lc_var->getCellID0() );
    pl_var.setCellID1( lc_var->getCellID1() );
    pl_var.setEDep( lc_var->getEDep() );
    pl_var.setTime( lc_var->getTime() );
    pl_var.setPathLength( lc_var->getPathLength() );
    pl_var.setQuality( lc_var->getQuality() );
    pl_var.setPosition( lc_var->getPosition() );
    pl_var.setMomentum( lc_var->getMomentum() );
    pl_var.setOverlay( lc_var->isOverlay() );
    pl_var.setProducedBySecondary( lc_var->isProducedBySecondary() );

    // Looping the LCIO::MCParticleCollection to pick Particle for Hits;
    CollectionsVec vec_mcp;
    vec_mcp = map_cols["MCParticle"];
    EVENT::LCCollection* hitcol_lc = (EVENT::LCCollection*) vec_mcp[0].first;
    plcio::MCParticleCollection* hitcol_pl = (plcio::MCParticleCollection*) vec_mcp[0].second;

    int index = -1;
    // search corresponding MCParticleCollection*;
    EVENT::MCParticle* mcptr_reg = lc_var->getMCParticle();
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
      pl_var.setMCParticle( hitcol_pl->at(index) );

  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_SimCalorimeterHit(EVENT::LCCollection* lc_col){
  plcio::SimCalorimeterHitCollection* pl_col = new plcio::SimCalorimeterHitCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::SimCalorimeterHit* lc_var = (EVENT::SimCalorimeterHit*) lc_col->getElementAt(i);
    plcio::SimCalorimeterHit pl_var = (plcio::SimCalorimeterHit) pl_col->create();

    pl_var.setCellID0( lc_var->getCellID0() );
    pl_var.setCellID1( lc_var->getCellID1() );
    pl_var.setEnergy( lc_var->getEnergy() );

    float plcio_v[3];
    const float* lcio_v = lc_var->getPosition();
    plcio_v[0] = lcio_v[0];
    plcio_v[1] = lcio_v[1];
    plcio_v[2] = lcio_v[2];
    pl_var.setPosition( plcio::FloatThree( plcio_v ) );

    // converting from lc_var to pl_var on the contribution variables;
    for( unsigned j=0, N=lc_var->getNMCContributions(); j<N; j++){
      plcio_v[0] = lc_var->getStepPosition(j)[0];
      plcio_v[1] = lc_var->getStepPosition(j)[1];
      plcio_v[2] = lc_var->getStepPosition(j)[2];

      plcio::ConstCaloHitContribution tmp(
	lc_var->getPDGCont(j), lc_var->getEnergyCont(j),
	lc_var->getTimeCont(j), plcio::FloatThree(plcio_v) 
      );
      pl_var.addContribution( tmp );
    }
  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_Cluster(EVENT::LCCollection* lc_col){
  plcio::ClusterCollection* pl_col = new plcio::ClusterCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::Cluster* lc_var = (EVENT::Cluster*) lc_col->getElementAt(i);
    plcio::Cluster pl_var = (plcio::Cluster) pl_col->create();

    pl_var.setType(lc_var->getType());
    pl_var.setEnergy(lc_var->getEnergy());
    pl_var.setEnergyError(lc_var->getEnergyError());
    pl_var.setPhi(lc_var->getIPhi());
    pl_var.setITheta(lc_var->getITheta());

    float plcio_v[3];
    const float* lcio_v = lc_var->getPosition();
    plcio_v[0] = lcio_v[0];
    plcio_v[1] = lcio_v[1];
    plcio_v[2] = lcio_v[2];
    pl_var.setPosition( plcio::FloatThree(plcio_v) );

    std::vector<float> lcio_vecfloat = lc_var->getPositionError();
    std::array<float, 6> plcio_arr6;
    for(unsigned j=0; j<6; j++){
      plcio_arr6[j] = lcio_vecfloat[j];
    }
    pl_var.setPositionError(plcio_arr6);

//    pl_var.(lc_var->getDirectionError());
//    pl_var.(lc_var->getShape());
//    pl_var.(lc_var->getParticleIDs());
//    pl_var.(lc_var->getClusters());
//    pl_var.(lc_var->getCalorimeterHits());
//    pl_var.(lc_var->getHitContributions());
//    pl_var.(lc_var->getSubdetectorEnergies());
  }

  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_Vertex(EVENT::LCCollection* lc_col){
  plcio::VertexCollection* pl_col = new plcio::VertexCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::Vertex* lc_var = (EVENT::Vertex*) lc_col->getElementAt(i);
    plcio::Vertex pl_var = (plcio::Vertex) pl_col->create();
  
    // Quest: data.primary is an int value, set by boolean number;
    pl_var.setPrimary( lc_var->isPrimary() );
    pl_var.setChi2( lc_var->getChi2() );
    pl_var.setProbability( lc_var->getProbability() );

    float plcio_v[3];
    const float* lcio_v = lc_var->getPosition();
    plcio_v[0] = lcio_v[0];
    plcio_v[1] = lcio_v[1];
    plcio_v[2] = lcio_v[2];
    pl_var.setPosition( plcio::FloatThree(plcio_v) );

    std::vector<float> vec_pra = lc_var->getParameters();
    for( unsigned j=0,M=vec_pra.size(); j<M; j++){
      pl_var.addParameter(vec_pra[j]);
    }

    // convert string into int(type code);
//    pl_var.setAlgorithmType( lc_var->getAlgorithmType() );
//    pl_var.setCovMatrix( lc_var->() );
//    pl_var.setAssociatedParticle( lc_var->() );
  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_TPCHit(EVENT::LCCollection* lc_col){
  plcio::TPCHitCollection* pl_col = new plcio::TPCHitCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::TPCHit* lc_var = (EVENT::TPCHit*) lc_col->getElementAt(i);
    plcio::TPCHit pl_var = (plcio::TPCHit) pl_col->create();

    pl_var.setCellID(lc_var->getCellID());
    pl_var.setTime(lc_var->getTime());
    pl_var.setCharge(lc_var->getCharge());
    pl_var.setQuality(lc_var->getQuality());

    for( unsigned j=0,M=lc_var->getNRawDataWords(); j<M; j++){
      pl_var.addRawDataWord(lc_var->getRawDataWord(j));
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
