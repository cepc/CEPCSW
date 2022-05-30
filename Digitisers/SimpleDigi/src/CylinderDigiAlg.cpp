#include "CylinderDigiAlg.h"

#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Vector3D.h"

#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

DECLARE_COMPONENT( CylinderDigiAlg )

CylinderDigiAlg::CylinderDigiAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc){
  // Input collections
  declareProperty("InputSimTrackerHitCollection", m_inputColHdls, "Handle of the Input SimTrackerHit collection");

  // Output collections
  declareProperty("OutputTrackerHitCollection", m_outputColHdls, "Handle of the output TrackerHit collection");
  declareProperty("TrackerHitAssociationCollection", m_assColHdls, "Handle of the Association collection between SimTrackerHit and TrackerHit");
}

StatusCode CylinderDigiAlg::initialize(){
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if(!m_geosvc){
    error() << "Failed to get the GeomSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  auto detector = m_geosvc->lcdd();
  if(!detector){
    error() << "Failed to get the Detector from GeomSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  std::string name = m_inputColHdls.objKey(); 
  debug() << "Readout name: " << name << endmsg;
  m_decoder = m_geosvc->getDecoder(name);
  if(!m_decoder){
    error() << "Failed to get the decoder. " << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "CylinderDigiAlg::initialized" << endmsg;
  return GaudiAlgorithm::initialize();
}


StatusCode CylinderDigiAlg::execute(){
  auto trkhitVec = m_outputColHdls.createAndPut();
  auto assVec = m_assColHdls.createAndPut();

  const edm4hep::SimTrackerHitCollection* STHCol = nullptr;
  try {
    STHCol = m_inputColHdls.get();
  }
  catch(GaudiException &e){
    debug() << "Collection " << m_inputColHdls.fullKey() << " is unavailable in event " << m_nEvt << endmsg;
    return StatusCode::SUCCESS;
  }
  if(STHCol->size()==0) return StatusCode::SUCCESS;
  debug() << m_inputColHdls.fullKey() << " has SimTrackerHit "<< STHCol->size() << endmsg;
  
  for(auto simhit : *STHCol){
    auto particle = simhit.getMCParticle();
    if(!particle.isAvailable()) continue;
    
    auto& mom0 = particle.getMomentum();
    double pt = sqrt(mom0.x*mom0.x+mom0.y*mom0.y);
    if(particle.isCreatedInSimulation()&&pt<0.01&&particle.isStopped()) continue;

    auto cellId = simhit.getCellID();
    int system  = m_decoder->get(cellId, "system");
    int chamber = m_decoder->get(cellId, "chamber");
    int layer   = m_decoder->get(cellId, "layer"  );
    auto& pos   = simhit.getPosition();
    auto& mom   = simhit.getMomentum();
    
    double phi = atan2(pos[1], pos[0]);
    double r   = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    double dphi = m_resRPhi/r;
    phi += randSvc()->generator(Rndm::Gauss(0, dphi))->shoot();
    double smearedX = r*cos(phi);
    double smearedY = r*sin(phi);
    double smearedZ = pos[2] + randSvc()->generator(Rndm::Gauss(0, m_resZ))->shoot();

    auto trkHit = trkhitVec->create();
    trkHit.setCellID(cellId);
    trkHit.setTime(simhit.getTime());
    trkHit.setEDep(simhit.getEDep());
    trkHit.setEdx (simhit.getEDep()/simhit.getPathLength());
    trkHit.setPosition (edm4hep::Vector3d(smearedX, smearedY, smearedZ));
    trkHit.setCovMatrix(std::array<float, 6>{m_resRPhi*m_resRPhi/2, 0, m_resRPhi*m_resRPhi/2, 0, 0, m_resZ*m_resZ});
    //trkHit.setType(CEPC::CYLINDER);
    trkHit.addToRawHits(simhit.getObjectID());
    debug() << "Hit " << simhit.id() << ": " << pos << " -> " << trkHit.getPosition() << "s:" << system << " c:" << chamber << " l:" << layer
	    << " pt = " << pt << " " << mom.x << " " << mom.y << " " << mom.z << endmsg;

    auto ass = assVec->create();

    float weight = 1.0;

    debug() <<" Set relation between " << " sim hit " << simhit.id() << " to tracker hit " << trkHit.id() << " with a weight of " << weight << endmsg;
    ass.setSim(simhit);
    ass.setRec(trkHit);
    ass.setWeight(weight);
  }

  m_nEvt++;

  return StatusCode::SUCCESS;
}

StatusCode CylinderDigiAlg::finalize(){
  info() << "Processed " << m_nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}
