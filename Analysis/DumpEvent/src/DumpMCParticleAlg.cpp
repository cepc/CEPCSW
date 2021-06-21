#include "DumpMCParticleAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include <math.h>

DECLARE_COMPONENT( DumpMCParticleAlg )

//------------------------------------------------------------------------------
DumpMCParticleAlg::DumpMCParticleAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
  declareProperty("MCParticleCollection", _inMCColHdl, "Handle of the Input MCParticle collection");

  m_thisName = name;
}

//------------------------------------------------------------------------------
StatusCode DumpMCParticleAlg::initialize(){
  info() << "Booking Ntuple" << endmsg;

  NTuplePtr nt1(ntupleSvc(), "MyTuples/MC");
  if ( !nt1 ) {
    m_tuple = ntupleSvc()->book("MyTuples/MC",CLID_ColumnWiseTuple,"MC truth");
    if ( 0 != m_tuple ) {
      m_tuple->addItem ("nmc",        m_nParticles, 0, 1000 ).ignore();
      m_tuple->addIndexedItem ("pdg",        m_nParticles, m_pdgID ).ignore();
      m_tuple->addIndexedItem ("genStatus",  m_nParticles, m_genStatus ).ignore();
      m_tuple->addIndexedItem ("simStatus",  m_nParticles, m_simStatus ).ignore();
      m_tuple->addIndexedItem ("charge",     m_nParticles, m_charge ).ignore();
      m_tuple->addIndexedItem ("time",       m_nParticles, m_time ).ignore();
      m_tuple->addIndexedItem ("mass",       m_nParticles, m_mass ).ignore();
      m_tuple->addIndexedItem ("vx",         m_nParticles, m_vx ).ignore();
      m_tuple->addIndexedItem ("vy",         m_nParticles, m_vy ).ignore();
      m_tuple->addIndexedItem ("vz",         m_nParticles, m_vz ).ignore();
      m_tuple->addIndexedItem ("px",         m_nParticles, m_px ).ignore();
      m_tuple->addIndexedItem ("py",         m_nParticles, m_py ).ignore();
      m_tuple->addIndexedItem ("pz",         m_nParticles, m_pz ).ignore();
      m_tuple->addIndexedItem ("d0",         m_nParticles, m_d0 ).ignore();
      m_tuple->addIndexedItem ("phi0",       m_nParticles, m_phi0 ).ignore();
      m_tuple->addIndexedItem ("omega",      m_nParticles, m_omega ).ignore();
      m_tuple->addIndexedItem ("z0",         m_nParticles, m_z0 ).ignore();
      m_tuple->addIndexedItem ("tanLambda",  m_nParticles, m_tanLambda ).ignore();
    }
    else { // did not manage to book the N tuple....
      fatal() << "Cannot bool MyTuples/MC " << endmsg;
      return StatusCode::FAILURE;
    }
  }
  else{
    m_tuple = nt1;
  }

  auto geomSvc = service<IGeomSvc>("GeomSvc");
  if(geomSvc){
    const dd4hep::Direction& field = geomSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
    m_field = field.z()/dd4hep::tesla;
    info() << "Magnetic field will obtain from GeomSvc = " << m_field << " tesla" << endmsg; 
  }
  else{
    info() << "Failed to find GeomSvc ..." << endmsg;
    info() << "Magnetic field will use what input through python option for this algorithm namse as Field, now " << m_field << " tesla" << endmsg;
  }

  _nEvt = 0;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode DumpMCParticleAlg::execute(){
  const edm4hep::MCParticleCollection* mcCols = nullptr;
  try {
    mcCols = _inMCColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inMCColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  if(mcCols){
    m_nParticles = 0;
    for(auto particle : *mcCols){
      m_pdgID[m_nParticles] = particle.getPDG();
      m_genStatus[m_nParticles] = particle.getGeneratorStatus();
      m_simStatus[m_nParticles] = particle.getSimulatorStatus();
      m_charge[m_nParticles] = particle.getCharge();
      m_time[m_nParticles] = particle.getTime();
      m_mass[m_nParticles] = particle.getMass();
      const edm4hep::Vector3d& vertex = particle.getVertex();
      m_vx[m_nParticles] = vertex.x;
      m_vy[m_nParticles] = vertex.y;
      m_vz[m_nParticles] = vertex.z;
      const edm4hep::Vector3f& momentum = particle.getMomentum();
      m_px[m_nParticles] = momentum.x;
      m_py[m_nParticles] = momentum.y;
      m_pz[m_nParticles] = momentum.z;

      HelixClass helix;
      float posV[3] = {vertex.x,vertex.y,vertex.z};
      float momV[3] = {momentum.x,momentum.y,momentum.z};
      helix.Initialize_VP(posV,momV,particle.getCharge(),m_field);
      float phiMC = helix.getPhi0();
      if(phiMC>CLHEP::pi) phiMC = phiMC - CLHEP::twopi;
      m_phi0[m_nParticles] = phiMC;
      m_d0[m_nParticles] = helix.getD0();
      m_omega[m_nParticles] = helix.getOmega();
      m_z0[m_nParticles] = helix.getZ0();
      m_tanLambda[m_nParticles] = helix.getTanLambda();
      m_nParticles++;
    }
    debug() << "MCParticle: " << m_nParticles <<endmsg;
  }
  m_tuple->write();
  _nEvt++;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode DumpMCParticleAlg::finalize(){
  debug() << "Finalizing..." << endmsg;

  return StatusCode::SUCCESS;
}
