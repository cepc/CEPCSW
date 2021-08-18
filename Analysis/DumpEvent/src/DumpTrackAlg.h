#ifndef DumpTrackAlg_h
#define DumpTrackAlg_h 1

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/Algorithm.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"

#include "GaudiKernel/NTuple.h"

class DumpTrackAlg : public Algorithm {
 public:
  // Constructor of this form must be provided
  DumpTrackAlg( const std::string& name, ISvcLocator* pSvcLocator );

  // Three mandatory member functions of any algorithm
  StatusCode initialize() override;
  StatusCode execute() override;
  StatusCode finalize() override;

 private:
  DataHandle<edm4hep::MCParticleCollection> _inMCColHdl{"MCParticle", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackCollection> _inTrackColHdl{"SiTracks", Gaudi::DataHandle::Reader, this};

  Gaudi::Property<double> m_field{this, "Field", 3.0};

  NTuple::Tuple*       m_tuple;
  NTuple::Item<long>   m_nTracks;
  NTuple::Array<float> m_x;
  NTuple::Array<float> m_y;
  NTuple::Array<float> m_z;
  NTuple::Array<float> m_px;
  NTuple::Array<float> m_py;
  NTuple::Array<float> m_pz;
  NTuple::Array<float> m_d0;
  NTuple::Array<float> m_phi0;
  NTuple::Array<float> m_omega;
  NTuple::Array<float> m_z0;
  NTuple::Array<float> m_tanLambda;
  NTuple::Array<float> m_sigma_d0;
  NTuple::Array<float> m_sigma_phi0;
  NTuple::Array<float> m_sigma_omega;
  NTuple::Array<float> m_sigma_z0;
  NTuple::Array<float> m_sigma_tanLambda;
  NTuple::Array<int>   m_nHitsVXD;
  NTuple::Array<int>   m_nHitsFTD;
  NTuple::Array<int>   m_nHitsSIT;
  NTuple::Array<int>   m_nHitsGAS;
  NTuple::Array<int>   m_nHitsSET;
  
  int _nEvt;
  std::string m_thisName;
};

#endif
