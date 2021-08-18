#ifndef DumpMCParticleAlg_h
#define DumpMCParticleAlg_h 1

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"

#include "GaudiKernel/NTuple.h"

class DumpMCParticleAlg : public Algorithm {
 public:
  // Constructor of this form must be provided
  DumpMCParticleAlg( const std::string& name, ISvcLocator* pSvcLocator );

  // Three mandatory member functions of any algorithm
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;

 private:
  DataHandle<edm4hep::MCParticleCollection> _inMCColHdl{"MCParticle", Gaudi::DataHandle::Reader, this};

  Gaudi::Property<double> m_field{this, "Field", 3.0};

  NTuple::Tuple*        m_tuple;
  NTuple::Item<long>     m_nParticles;
  NTuple::Array<int>    m_pdgID;
  NTuple::Array<int>    m_genStatus;
  NTuple::Array<int>    m_simStatus;
  NTuple::Array<float>  m_charge;
  NTuple::Array<float>  m_time;
  NTuple::Array<double> m_mass;
  NTuple::Array<double> m_vx;
  NTuple::Array<double> m_vy;
  NTuple::Array<double> m_vz;
  NTuple::Array<float>  m_px;
  NTuple::Array<float>  m_py;
  NTuple::Array<float>  m_pz;
  NTuple::Array<float>  m_d0;
  NTuple::Array<float>  m_phi0;
  NTuple::Array<float>  m_omega;
  NTuple::Array<float>  m_z0;
  NTuple::Array<float>  m_tanLambda;
  
  int _nEvt;
  std::string m_thisName;
};

#endif
