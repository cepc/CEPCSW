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
std::string LCIO2Plcio::CollName;

void lp_info(std::string s){
  printf("[LCIO2Plcio]:\t\t%s.\n", &s[0]);
}

template<typename T>
static plcio::FloatThree FloatThreeFROMConstPtr(const T* vpos){
  float tmp[3];
  for(unsigned i=0; i<3; i++){
    tmp[i] = vpos[i];
  }
  return plcio::FloatThree(tmp);
}

template<typename T>
static plcio::FloatThree FloatThreeFROMFloatVec(std::vector<T> vec){
  float tmp[3];
  for(unsigned i=0; i<3; i++){
    tmp[i] = vec[i];
  }
  return plcio::FloatThree(tmp);
}

template<typename T>
static plcio::DoubleThree DoubleThreeFROMConstPtr(const T* vpos){
  double tmp[3];
  for(unsigned i=0; i<3; i++)
    tmp[i] = vpos[i];
  return plcio::DoubleThree(tmp);
}

std::array<float, 6> vec6_2_arr6(std::vector<float> vec){
  std::array<float, 6> arr;
  for(unsigned i=0; i<6; i++){
    arr[i] = vec[i];
  }
  return arr;
}

LCIO2Plcio::LCIO2Plcio(){
  map_cvt.insert(std::make_pair<std::string, fptr>("MCParticle", Convertor_MCParticle));
  map_cvt.insert(std::make_pair<std::string, fptr>("LCRunHeader", Convertor_LCRunHeader));
  map_cvt.insert(std::make_pair<std::string, fptr>("SimTrackerHit", Convertor_SimTrackerHit));
  map_cvt.insert(std::make_pair<std::string, fptr>("SimCalorimeterHit", Convertor_SimCalorimeterHit));
  map_cvt.insert(std::make_pair<std::string, fptr>("Cluster", Convertor_Cluster));
  map_cvt.insert(std::make_pair<std::string, fptr>("Track", Convertor_Track));
  map_cvt.insert(std::make_pair<std::string, fptr>("TrackerHit", Convertor_TrackerHit));
  map_cvt.insert(std::make_pair<std::string, fptr>("TPCHit", Convertor_TPCHit));
  map_cvt.insert(std::make_pair<std::string, fptr>("ReconstructedParticle", Convertor_ReconstructedParticle));
  map_cvt.insert(std::make_pair<std::string, fptr>("ParticleID", Convertor_ParticleID));
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
  pl_var.setMomentum( FloatThreeFROMConstPtr(lc_var->getMomentum()) );
  pl_var.setMomentumAtEndpoint( FloatThreeFROMConstPtr(lc_var->getMomentumAtEndpoint()) );
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
    pl_var.setPosition( FloatThreeFROMConstPtr(lc_var->getPosition()) );

    // converting from lc_var to pl_var on the contribution variables;
    for( unsigned j=0, N=lc_var->getNMCContributions(); j<N; j++){
      plcio::ConstCaloHitContribution tmp(
	lc_var->getPDGCont(j), lc_var->getEnergyCont(j),
	lc_var->getTimeCont(j), FloatThreeFROMConstPtr(lc_var->getStepPosition(j)) 
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
    pl_var.setPosition( FloatThreeFROMConstPtr(lc_var->getPosition()) );
    pl_var.setPositionError( vec6_2_arr6(lc_var->getPositionError()) );
    pl_var.setDirectionError( FloatThreeFROMFloatVec(lc_var->getDirectionError()) );

    std::vector<float> lcio_seq0 = lc_var->getShape();
    for(unsigned j=0; j<lcio_seq0.size(); j++){
      pl_var.addShap( lcio_seq0[j] );
    }

    lcio_seq0 = lc_var->getHitContributions();
    for(unsigned j=0; j<lcio_seq0.size(); j++){
      pl_var.addHitContribution( lcio_seq0[j] );
    }

    lcio_seq0 = lc_var->getSubdetectorEnergies();
    for(unsigned j=0; j<lcio_seq0.size(); j++){
      pl_var.addSubdetectorEnergie( lcio_seq0[j] );
    }

    EVENT::ParticleIDVec lcio_seq1 = lc_var->getParticleIDs();
    for(unsigned j=0; j<lcio_seq1.size(); j++){
      EVENT::ParticleID* lc_locx = lcio_seq1[j];

      pl_var.addParticleID( plcio::ConstParticleID(
        lc_locx->getType(), 
	lc_locx->getPDG(),
	lc_locx->getAlgorithmType(), 
	lc_locx->getLikelihood()
      ) );
    }

    EVENT::CalorimeterHitVec lcio_seq2 = lc_var->getCalorimeterHits();
    for(unsigned j=0; j<lcio_seq2.size(); j++){
      EVENT::CalorimeterHit* lc_locx = lcio_seq2[j];

      podio::CollectionIDTable col2id;
      plcio::ObjectID tmp;
      tmp.index = i;
      tmp.collectionID = col2id.collectionID(CollName);

      pl_var.addHit( plcio::ConstCalorimeterHit(
        lc_locx->getCellID0(),
        lc_locx->getCellID1(),
        lc_locx->getEnergy(),
        lc_locx->getEnergyError(),
        lc_locx->getTime(),
        FloatThreeFROMConstPtr(lc_locx->getPosition()),
        lc_locx->getType(),
        tmp
      ) );
    }
  }

  for(unsigned i=0,N=lc_col->getNumberOfElements(); i<N; i++){
    EVENT::Cluster* lc_var = (EVENT::Cluster*) lc_col->getElementAt(i);
    plcio::Cluster pl_var = pl_col->at(i);
    EVENT::ClusterVec lcio_seq3 = lc_var->getClusters();

    for(unsigned j=0; j<lcio_seq3.size(); j++){
      for(unsigned k=0,K=lc_col->getNumberOfElements(); k<K; k++){
        if(lcio_seq3[j] == lc_col->getElementAt(k)){
          pl_var.addCluster(pl_col->at(k));
	  break;
	}
      }
    }
  }

  return pl_col;
}

// QUEST::isAvailable dows not analyze;
podio::CollectionBase* LCIO2Plcio::Convertor_ParticleID(EVENT::LCCollection* lc_col){
  plcio::ParticleIDCollection* pl_col = new plcio::ParticleIDCollection();

  for( unsigned i=0,N=lc_col->getNumberOfElements() ; i< N ; ++i){
    EVENT::ParticleID* lc_var = (EVENT::ParticleID*) lc_col->getElementAt(i);
    plcio::ParticleID pl_var = (plcio::ParticleID) pl_col->create();

    pl_var.setType( lc_var->getType() );
    pl_var.setPDG( lc_var->getPDG() );
    pl_var.setAlgorythmType( lc_var->getAlgorithmType() );
    pl_var.setLikelihood( lc_var->getLikelihood() );
    EVENT::FloatVec vec_float = lc_var->getParameters();
    for(unsigned j=0; j<vec_float.size(); j++){
      pl_var.addParameter(vec_float[j]);
    }
  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_ReconstructedParticle(EVENT::LCCollection* lc_col){
  plcio::ReconstructedParticleCollection* pl_col = new plcio::ReconstructedParticleCollection();

  for(unsigned i=0,N=lc_col->getNumberOfElements(); i<N; i++){
    EVENT::ReconstructedParticle* lc_var = (EVENT::ReconstructedParticle*) lc_col->getElementAt(i);
    plcio::ReconstructedParticle pl_var = (plcio::ReconstructedParticle) pl_col->create();

    pl_var.setType( lc_var->getType() );
    pl_var.setEnergy( lc_var->getEnergy() );
    pl_var.setCharge( lc_var->getCharge() );
    pl_var.setMass( lc_var->getMass() );
    pl_var.setGoodnessOfPID( lc_var->getGoodnessOfPID() );

    pl_var.setMomentum( FloatThreeFROMConstPtr(lc_var->getMomentum()) );
    pl_var.setReferencePoint( FloatThreeFROMConstPtr(lc_var->getReferencePoint()) );

    std::vector<float> vec = lc_var->getCovMatrix();
    for(unsigned j=0,N=vec.size(); j<N; j++){
      pl_var.setCovMatrix( j, vec[j] );
    }

    // QUEST: boolean to int: isPrimary, str2int: getAlgorithmType;
    // pl_var.setStartVertex( lc_var->getStartVertex() ) 
    // ignorant;
    std::array<float, 6> arr_6;
    EVENT::FloatVec fvec = lc_var->getStartVertex()->getCovMatrix();
    for(unsigned j=0; j<6; j++){
      arr_6[j] = fvec[j];
    }
    pl_var.setStartVertex( plcio::ConstVertex(
      lc_var->getStartVertex()->isPrimary(),
      lc_var->getStartVertex()->getChi2(),
      lc_var->getStartVertex()->getProbability(),
      FloatThreeFROMConstPtr( lc_var->getStartVertex()->getPosition() ),
      arr_6,
      0
      //lc_var->getStartVertex()->getAlgorithmType()
    ) );

    //pl_var.setParticleIDUsed( lc_var->getParticleIDUsed() );
    EVENT::ParticleID* lc_locx = lc_var->getParticleIDUsed();
    CollectionsVec vec_cols = map_cols["ParticleID"];
    EVENT::LCCollection* lc_mapcol = (EVENT::LCCollection*)vec_cols[0].first;
    plcio::ParticleIDCollection* pl_mapcol = (plcio::ParticleIDCollection*)vec_cols[0].second;

    for(unsigned k=0, LCsize=lc_mapcol->getNumberOfElements(); k<LCsize; k++){
      if(lc_locx == lc_mapcol->getElementAt(k))
        pl_var.setParticleIDUsed(pl_mapcol->at(k));
    }
    

    //pl_var.addCluster();
    std::vector<EVENT::Cluster*> vec_clust = lc_var->getClusters();
    for(unsigned j=0; j<vec_clust.size(); j++){
      CollectionsVec vec_cols = map_cols["Cluster"];
      EVENT::LCCollection* lc_mapcol = (EVENT::LCCollection*)vec_cols[0].first;
      plcio::ClusterCollection* pl_mapcol = (plcio::ClusterCollection*)vec_cols[0].second;

      for(unsigned k=0; k<lc_mapcol->getNumberOfElements(); k++){
        if( vec_clust[j] == lc_mapcol->getElementAt(k) )
          pl_var.addCluster( pl_mapcol->at(k) );
      }
    }

    //pl_var.addTrack();
    EVENT::TrackVec vec_track = lc_var->getTracks();
    for(unsigned j=0; j<vec_track.size(); j++){
      CollectionsVec vec_cols = map_cols["Track"];
      EVENT::LCCollection* lc_mapcol = (EVENT::LCCollection*)vec_cols[0].first;
      plcio::TrackCollection* pl_mapcol = (plcio::TrackCollection*)vec_cols[0].second;

      for(unsigned k=0; k<lc_mapcol->getNumberOfElements(); k++)
        if(vec_track[j] == lc_mapcol->getElementAt(k))
          pl_var.addTrack(pl_mapcol->at(k));
    }

    //pl_var.addParticleID( lc->getParticleIDs);
    EVENT::ParticleIDVec vec_ParticleID = lc_var->getParticleIDs();
    for(unsigned j=0; j<vec_ParticleID.size(); j++){
      CollectionsVec vec_cols = map_cols["ParticleID"];
      EVENT::LCCollection* lc_mapcol = (EVENT::LCCollection*)vec_cols[0].first;
      plcio::ParticleIDCollection* pl_mapcol = (plcio::ParticleIDCollection*)vec_cols[0].second;

      for(unsigned k=0,M=lc_mapcol->getNumberOfElements(); k<M; k++){
        if(vec_ParticleID[j] == lc_mapcol->getElementAt(k))
          pl_var.addParticleID( pl_mapcol->at(k) );
      }
    }
  }
  //pl_var.addParticle();
  for(unsigned i=0, N=lc_col->getNumberOfElements(); i<N; i++){
    EVENT::ReconstructedParticle* lc_var = (EVENT::ReconstructedParticle*)lc_col->getElementAt(i);
    EVENT::ReconstructedParticleVec vec_RecPtc = lc_var->getParticles();

    for(unsigned j=0; j<vec_RecPtc.size(); j++){
      for(unsigned k=0, M=lc_col->getNumberOfElements(); k<M; k++){
        if(vec_RecPtc[j] == lc_col->getElementAt(k))
          pl_col->at(i).addParticle( pl_col->at(k) );
      }
    }
  }

  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_TrackerHit(EVENT::LCCollection* lc_col){
  plcio::TrackerHitCollection* pl_col = new plcio::TrackerHitCollection();
  
  for(unsigned i=0,N=lc_col->getNumberOfElements(); i<N; i++){
    EVENT::TrackerHit* lc_var = (EVENT::TrackerHit*) lc_col->getElementAt(i);
    plcio::TrackerHit  pl_var = (plcio::TrackerHit) pl_col->create();
    
    pl_var.setCellID0(lc_var->getCellID0());
    pl_var.setCellID1(lc_var->getCellID1());

    pl_var.setType(lc_var->getType());
    pl_var.setQuality(lc_var->getQuality());
    pl_var.setTime(lc_var->getTime());
    pl_var.setEDep(lc_var->getEDep());
    pl_var.setEDepError(lc_var->getEDepError());
    pl_var.setEdx(lc_var->getdEdx());
    pl_var.setPosition( DoubleThreeFROMConstPtr(lc_var->getPosition()));
    pl_var.setCovMatrix( vec6_2_arr6(lc_var->getCovMatrix()));
  }
  return pl_col;
}

podio::CollectionBase* LCIO2Plcio::Convertor_Track(EVENT::LCCollection* lc_col){
  plcio::TrackCollection* pl_col = new plcio::TrackCollection();

  for(unsigned i=0,N=lc_col->getNumberOfElements(); i<N; i++){
    EVENT::Track* lc_var = (EVENT::Track*) lc_col->getElementAt(i);
    plcio::Track pl_var = (plcio::Track) pl_col->create();

    pl_var.setType( lc_var->getType() );
    pl_var.setChi2( lc_var->getChi2() );
    pl_var.setNdf( lc_var->getNdf() );
    pl_var.setDEdx( lc_var->getdEdx() );
    pl_var.setDEdxError( lc_var->getdEdxError() );
    pl_var.setRadiusOfInnermostHit( lc_var->getRadiusOfInnermostHit() );

    //pl_var.addTrackerHit( lc_var->getTrackerHits() );
    //rely on TrackerHits collection;
    //corresoponding with TrackerHit collection; 
    EVENT::TrackerHitVec lcio_seq1 = lc_var->getTrackerHits();
    for(unsigned j=0; j<lcio_seq1.size(); j++){
      EVENT::TrackerHit* lc_locx = lcio_seq1[j];
      CollectionsVec vec_cols = map_cols["TrackerHit"];
      EVENT::LCCollection* lc_mapcol = (EVENT::LCCollection*)vec_cols[0].first;
      plcio::TrackerHitCollection* pl_mapcol = (plcio::TrackerHitCollection*) vec_cols[0].second;

      for(unsigned k=0; k<lc_mapcol->getNumberOfElements(); k++){
        if( lc_locx == lc_mapcol->getElementAt(k) )
          pl_var.addTrackerHit( pl_mapcol->at(k) );
      } 
    }

    //pl_var.( lc_var->getSubdetectorHitNumbers() );
    std::vector<int> lcio_IntVec1 = lc_var->getSubdetectorHitNumbers();
    for(unsigned j=0; j<lcio_IntVec1.size(); j++){
      pl_var.addSubDetectorHitNumber(lcio_IntVec1[j]);
    }

    //pl_var.( lc_var->getTrackStates() );
    EVENT::TrackStateVec lcio_seq3 = lc_var->getTrackStates();
    for(unsigned j=0; j<lcio_seq3.size(); j++){
      EVENT::TrackState* lc_locx = lcio_seq3[j];
      std::array<float, 15> tmp_covM;
      for(unsigned k=0; k<15; k++){
        tmp_covM[k] = lc_locx->getCovMatrix()[k];
      }

      plcio::TrackState tmp;
      tmp.D0 = lc_locx->getD0();
      tmp.Z0 = lc_locx->getZ0();
      tmp.covMatrix = tmp_covM; 
      tmp.location = lc_locx->getLocation();
      tmp.omega = lc_locx->getOmega();
      tmp.phi = lc_locx->getPhi();
      tmp.referencePoint = FloatThreeFROMConstPtr( lc_locx->getReferencePoint() );
      tmp.tanLambda = lc_locx->getTanLambda();
      pl_var.addTrackState( tmp );
    }
  }

  //pl_var.( lc_var->getTracks() );
  for(unsigned i=0; i<lc_col->getNumberOfElements(); i++){
    EVENT::Track* lc_var = (EVENT::Track*)lc_col->getElementAt(i);
    EVENT::TrackVec lcio_seq2 = lc_var->getTracks();

    for(unsigned j=0; j<lcio_seq2.size(); j++){
      EVENT::Track* lc_locx = lcio_seq2[j];
      for(unsigned k=0; k<lc_col->getNumberOfElements(); k++){
        if( lc_locx == lc_col->getElementAt(k) )
          pl_col->at(i).addTrack(pl_col->at(k));
      }
    }
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
