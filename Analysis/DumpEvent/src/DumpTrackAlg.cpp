#include "DumpTrackAlg.h"

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

DECLARE_COMPONENT( DumpTrackAlg )

//------------------------------------------------------------------------------
DumpTrackAlg::DumpTrackAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
  declareProperty("MCParticleCollection", _inMCColHdl, "Handle of the Input MCParticle collection");
  declareProperty("TrackCollection", _inTrackColHdl, "Handle of the Input Track collection from CEPCSW");

  m_thisName = name;
}

//------------------------------------------------------------------------------
StatusCode DumpTrackAlg::initialize(){
  info() << "Booking Ntuple" << endmsg;

  NTuplePtr nt1(ntupleSvc(), "MyTuples/Track"+m_thisName);
  if ( !nt1 ) {
    m_tuple = ntupleSvc()->book("MyTuples/Track"+m_thisName,CLID_ColumnWiseTuple,"Tracking result");
    if ( 0 != m_tuple ) {
      m_tuple->addItem ("ntrk",      m_nTracks, 0, 500 ).ignore();
      m_tuple->addIndexedItem ("x",         m_nTracks, m_x ).ignore();
      m_tuple->addIndexedItem ("y",         m_nTracks, m_y ).ignore();
      m_tuple->addIndexedItem ("z",         m_nTracks, m_z ).ignore();
      m_tuple->addIndexedItem ("px",        m_nTracks, m_px ).ignore();
      m_tuple->addIndexedItem ("py",        m_nTracks, m_py ).ignore();
      m_tuple->addIndexedItem ("pz",        m_nTracks, m_pz ).ignore();
      m_tuple->addIndexedItem ("d0",        m_nTracks, m_d0 ).ignore();
      m_tuple->addIndexedItem ("phi0",      m_nTracks, m_phi0 ).ignore();
      m_tuple->addIndexedItem ("omega",     m_nTracks, m_omega ).ignore();
      m_tuple->addIndexedItem ("z0",        m_nTracks, m_z0 ).ignore();
      m_tuple->addIndexedItem ("tanLambda", m_nTracks, m_tanLambda ).ignore();
      m_tuple->addIndexedItem ("sigma_d0",        m_nTracks, m_sigma_d0 ).ignore();
      m_tuple->addIndexedItem ("sigma_phi0",      m_nTracks, m_sigma_phi0 ).ignore();
      m_tuple->addIndexedItem ("sigma_omega",     m_nTracks, m_sigma_omega ).ignore();
      m_tuple->addIndexedItem ("sigma_z0",        m_nTracks, m_sigma_z0 ).ignore();
      m_tuple->addIndexedItem ("sigma_tanLambda", m_nTracks, m_sigma_tanLambda ).ignore();
      m_tuple->addIndexedItem ("nvxd", m_nTracks, m_nHitsVXD ).ignore();
      m_tuple->addIndexedItem ("nftd", m_nTracks, m_nHitsFTD ).ignore();
      m_tuple->addIndexedItem ("nsit", m_nTracks, m_nHitsSIT ).ignore();
      m_tuple->addIndexedItem ("ngas", m_nTracks, m_nHitsGAS ).ignore();
      m_tuple->addIndexedItem ("nset", m_nTracks, m_nHitsSET ).ignore();
    }
    else { // did not manage to book the N tuple....
      fatal() << "Cannot book MyTuples/Track"+m_thisName <<endmsg; 
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
StatusCode DumpTrackAlg::execute(){
  const edm4hep::TrackCollection* trackCols = nullptr;
  try {
    trackCols = _inTrackColHdl.get();
  }
  catch ( GaudiException &e ) {
    debug() << "Collection " << _inTrackColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
  }

  if(trackCols){
    m_nTracks = 0;
    for(auto track : *trackCols){
      // since possible more than one location=1 TrackState (not deleted in reconstruction), always use last one 
      for(std::vector<edm4hep::TrackState>::const_iterator it=track.trackStates_end()-1; it!=track.trackStates_begin()-1; it--){
	edm4hep::TrackState trackState = *it;
	if(trackState.location!=1)continue;
	m_d0[m_nTracks]              = trackState.D0;
	m_phi0[m_nTracks]            = trackState.phi;
	m_omega[m_nTracks]           = trackState.omega;
	m_z0[m_nTracks]              = trackState.Z0;
	m_tanLambda[m_nTracks]       = trackState.tanLambda;
	m_sigma_d0[m_nTracks]        = std::sqrt(trackState.covMatrix[0]);
	m_sigma_phi0[m_nTracks]      = std::sqrt(trackState.covMatrix[2]);
	m_sigma_omega[m_nTracks]     = std::sqrt(trackState.covMatrix[5]);
	m_sigma_z0[m_nTracks]        = std::sqrt(trackState.covMatrix[9]);
        m_sigma_tanLambda[m_nTracks] = std::sqrt(trackState.covMatrix[14]);
	HelixClass helix_fit;
        helix_fit.Initialize_Canonical(trackState.phi,trackState.D0,trackState.Z0,trackState.omega,trackState.tanLambda,m_field);
	m_x[m_nTracks] = helix_fit.getReferencePoint()[0];
	m_y[m_nTracks] = helix_fit.getReferencePoint()[1];
	m_z[m_nTracks] = helix_fit.getReferencePoint()[2];
        m_px[m_nTracks] = helix_fit.getMomentum()[0];
        m_py[m_nTracks] = helix_fit.getMomentum()[1];
        m_pz[m_nTracks] = helix_fit.getMomentum()[2];
	//info() << "size = " << track.subDetectorHitNumbers_size() << endmsg;
	//for(int ii=0;ii<track.subDetectorHitNumbers_size();ii++){
	//  std::cout << track.getSubDetectorHitNumbers(ii) << " ";
	//}
	//std::cout << std::endl;
	if(track.subDetectorHitNumbers_size()>=5){
	  m_nHitsVXD[m_nTracks] = track.getSubDetectorHitNumbers(0);
	  m_nHitsFTD[m_nTracks] = track.getSubDetectorHitNumbers(1);
	  m_nHitsSIT[m_nTracks] = track.getSubDetectorHitNumbers(2);
	  m_nHitsGAS[m_nTracks] = track.getSubDetectorHitNumbers(3);
	  m_nHitsSET[m_nTracks] = track.getSubDetectorHitNumbers(4);
	}
	else{
	  m_nHitsVXD[m_nTracks] = 0;
          m_nHitsSIT[m_nTracks] = 0;
          m_nHitsSET[m_nTracks] = 0;
          m_nHitsFTD[m_nTracks] = 0;
          m_nHitsGAS[m_nTracks] = 0;
	}
	m_nTracks++;
	break;
      }
    }
  }
  m_tuple->write();
  _nEvt++;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode DumpTrackAlg::finalize(){
  debug() << "Finalizing..." << endmsg;

  return StatusCode::SUCCESS;
}
