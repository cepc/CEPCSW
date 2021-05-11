/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "GearSvc/IGearSvc.h"
#include "PandoraMatrixAlg.h"
#include "EventSeeder/IEventSeeder.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CaloHitContribution.h"
#include "edm4hep/ClusterConst.h"
#include "UTIL/ILDConf.h"
#include <cmath>
#include <algorithm>
#include "gear/BField.h"
#include <gear/GEAR.h>

#include "LCContent.h"

#include "TInterpreter.h"

pandora::Pandora* PandoraMatrixAlg::m_pPandora=0;

DECLARE_COMPONENT( PandoraMatrixAlg )

template<typename T ,typename T1>
StatusCode getCol(T & t, T1 & t1)
{
    try {
        t1 = t.get();
    }
    catch ( GaudiException &e ) {
        std::cout << "Collection " << t.fullKey() << " is unavailable in event "  << std::endl;
    }
    return StatusCode::SUCCESS;
}



PandoraMatrixAlg::PandoraMatrixAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
 m_CollectionMaps = new CollectionMaps();
  
 declareProperty("WriteClusterCollection"              , m_ClusterCollection_w,               "Handle of the ClusterCollection               output collection" );
 declareProperty("WriteReconstructedParticleCollection", m_ReconstructedParticleCollection_w, "Handle of the ReconstructedParticleCollection output collection" );
 declareProperty("WriteVertexCollection"               , m_VertexCollection_w,                "Handle of the VertexCollection                output collection" );

}


void PandoraMatrixAlg::FinaliseSteeringParameters(ISvcLocator* svcloc)
{
    // ATTN: This function seems to be necessary for operations that cannot easily be performed at construction of the processor,
    // when the steering file is parsed e.g. the call to GEAR to get the inner bfield
    
    m_caloHitCreatorSettings.m_absorberRadLengthECal = m_geometryCreatorSettings.m_absorberRadLengthECal;
    m_caloHitCreatorSettings.m_absorberIntLengthECal = m_geometryCreatorSettings.m_absorberIntLengthECal;
    m_caloHitCreatorSettings.m_absorberRadLengthHCal = m_geometryCreatorSettings.m_absorberRadLengthHCal;
    m_caloHitCreatorSettings.m_absorberIntLengthHCal = m_geometryCreatorSettings.m_absorberIntLengthHCal;
    m_caloHitCreatorSettings.m_absorberRadLengthOther = m_geometryCreatorSettings.m_absorberRadLengthOther;
    m_caloHitCreatorSettings.m_absorberIntLengthOther = m_geometryCreatorSettings.m_absorberIntLengthOther;

    m_caloHitCreatorSettings.m_hCalEndCapInnerSymmetryOrder = m_geometryCreatorSettings.m_hCalEndCapInnerSymmetryOrder;
    m_caloHitCreatorSettings.m_hCalEndCapInnerPhiCoordinate = m_geometryCreatorSettings.m_hCalEndCapInnerPhiCoordinate;
    
    m_trackCreatorSettings.m_prongSplitVertexCollections = m_trackCreatorSettings.m_prongVertexCollections;
    m_trackCreatorSettings.m_prongSplitVertexCollections.insert(m_trackCreatorSettings.m_prongSplitVertexCollections.end(),
        m_trackCreatorSettings.m_splitVertexCollections.begin(), m_trackCreatorSettings.m_splitVertexCollections.end());
    
    //m_settings.m_innerBField = marlin::Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z();
    IGearSvc*  iSvc = 0;
    StatusCode sc = svcloc->service("GearSvc", iSvc, false);
    if ( !sc ) 
    {
        throw "Failed to find GearSvc ...";
    }
    gear::GearMgr* _GEAR = iSvc->getGearMgr();
    m_settings.m_innerBField = _GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z();
    std::cout<<"m_innerBField="<<m_settings.m_innerBField<<std::endl;    
    m_mcParticleCreatorSettings.m_bField = m_settings.m_innerBField;
}




StatusCode PandoraMatrixAlg::initialize()
{

  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  std::cout<<"hi init PandoraMatrixAlg"<<std::endl;

  std::string s_output =m_AnaOutput; 
  m_fout = new TFile(s_output.c_str(),"RECREATE"); 
  m_tree = new TTree("evt","tree");
  m_tree->Branch("m_pReco_PID"   , &m_pReco_PID);
  m_tree->Branch("m_pReco_mc_id" , &m_pReco_mc_id);
  m_tree->Branch("m_pReco_mc_weight", &m_pReco_mc_weight);
  m_tree->Branch("m_pReco_mass"  , &m_pReco_mass);
  m_tree->Branch("m_pReco_energy", &m_pReco_energy);
  m_tree->Branch("m_pReco_px"    , &m_pReco_px);
  m_tree->Branch("m_pReco_py"    , &m_pReco_py);
  m_tree->Branch("m_pReco_pz"    , &m_pReco_pz);
  m_tree->Branch("m_pReco_charge", &m_pReco_charge);

  m_tree->Branch("m_mc_id"    , &m_mc_id);
  m_tree->Branch("m_mc_p_size", &m_mc_p_size);
  m_tree->Branch("m_mc_pid"   , &m_mc_pid   );
  m_tree->Branch("m_mc_mass"  , &m_mc_mass  );
  m_tree->Branch("m_mc_px"    , &m_mc_px    );
  m_tree->Branch("m_mc_py"    , &m_mc_py    );
  m_tree->Branch("m_mc_pz"    , &m_mc_pz    );
  m_tree->Branch("m_mc_charge", &m_mc_charge);
  m_tree->Branch("m_hasConversion", &m_hasConversion);


  m_tree->Branch("m_hits_x", &m_hits_x);
  m_tree->Branch("m_hits_y", &m_hits_y);
  m_tree->Branch("m_hits_z", &m_hits_z);
  m_tree->Branch("m_hits_E", &m_hits_E);

  for ( const auto& col : m_readCols ) {
      auto seperater = col.find(':');
      std::string colType = col.substr(0, seperater);
      std::string colName = col.substr(seperater+1);
      m_collections[colName] = colType;

      if ( colType == "MCParticle" ) {
          m_dataHandles[colName] =
              new DataHandle<edm4hep::MCParticleCollection>(colName, Gaudi::DataHandle::Reader, this);
      }
      else if ( colType == "Track" ) {
          m_dataHandles[colName] =
              new DataHandle<edm4hep::TrackCollection>(colName, Gaudi::DataHandle::Reader, this);
      }
      else if ( colType == "CalorimeterHit" ) {
          m_dataHandles[colName] =
              new DataHandle<edm4hep::CalorimeterHitCollection>(colName, Gaudi::DataHandle::Reader, this);
      }
      else if ( colType == "Vertex" ) {
          m_dataHandles[colName] =
              new DataHandle<edm4hep::VertexCollection>(colName, Gaudi::DataHandle::Reader, this);
      }
      else if ( colType == "MCRecoTrackerAssociation" ) {
          m_dataHandles[colName] =
              new DataHandle<edm4hep::MCRecoTrackerAssociationCollection>(colName, Gaudi::DataHandle::Reader, this);
      }
      else if ( colType == "MCRecoCaloAssociation" ) {
          m_dataHandles[colName] =
              new DataHandle<edm4hep::MCRecoCaloAssociationCollection>(colName, Gaudi::DataHandle::Reader, this);
      }
      else {
            error() << "invalid collection type: " << colType << endmsg;
            return StatusCode::FAILURE;
      }
  }

  // XML file
  m_settings.m_pandoraSettingsXmlFile =  m_PandoraSettingsXmlFile ; 
  // Hadronic energy non-linearity correction
  m_settings.m_inputEnergyCorrectionPoints = m_InputEnergyCorrectionPoints;
  m_settings.m_outputEnergyCorrectionPoints = m_OutputEnergyCorrectionPoints;
  // B-field parameters
  m_settings.m_muonBarrelBField = m_MuonBarrelBField; 
  m_settings.m_muonEndCapBField = m_MuonEndCapBField;
  
  m_trackCreatorSettings.m_trackCollections = m_TrackCollections ; 
  m_trackCreatorSettings.m_kinkVertexCollections = m_KinkVertexCollections; 
  m_trackCreatorSettings.m_prongVertexCollections = m_ProngVertexCollections;
  m_trackCreatorSettings.m_splitVertexCollections = m_SplitVertexCollections;
  m_trackCreatorSettings.m_v0VertexCollections = m_V0VertexCollections; 
  
  m_caloHitCreatorSettings.m_eCalCaloHitCollections = m_ECalCaloHitCollections;
  m_caloHitCreatorSettings.m_hCalCaloHitCollections = m_HCalCaloHitCollections;
  m_caloHitCreatorSettings.m_lCalCaloHitCollections = m_LCalCaloHitCollections;
  m_caloHitCreatorSettings.m_lHCalCaloHitCollections = m_LHCalCaloHitCollections;
  m_caloHitCreatorSettings.m_muonCaloHitCollections = m_MuonCaloHitCollections; 
  m_mcParticleCreatorSettings.m_mcParticleCollections = m_MCParticleCollections;
  m_mcParticleCreatorSettings.m_CaloHitRelationCollections = m_RelCaloHitCollections; 
  m_mcParticleCreatorSettings.m_TrackRelationCollections = m_RelTrackCollections;
  
  
  // Absorber properties
  m_geometryCreatorSettings.m_absorberRadLengthECal = m_AbsorberRadLengthECal;
  m_geometryCreatorSettings.m_absorberIntLengthECal = m_AbsorberIntLengthECal;
  m_geometryCreatorSettings.m_absorberRadLengthHCal = m_AbsorberRadLengthHCal;
  m_geometryCreatorSettings.m_absorberIntLengthHCal = m_AbsorberIntLengthHCal;
  m_geometryCreatorSettings.m_absorberRadLengthOther = m_AbsorberRadLengthOther;
  m_geometryCreatorSettings.m_absorberIntLengthOther = m_AbsorberIntLengthOther;
  
  // Name of PFO collection written by GaudiPandora
  
  m_pfoCreatorSettings.m_clusterCollectionName = m_ClusterCollectionName;// not used  
  m_pfoCreatorSettings.m_pfoCollectionName = m_PFOCollectionName;//
  m_pfoCreatorSettings.m_startVertexCollectionName = m_StartVertexCollectionName; //
  m_pfoCreatorSettings.m_startVertexAlgName = m_StartVertexAlgorithmName;//
   
  m_pfoCreatorSettings.m_emStochasticTerm = m_EMStochasticTerm;
  m_pfoCreatorSettings.m_hadStochasticTerm = m_HadStochasticTerm;
  m_pfoCreatorSettings.m_emConstantTerm = m_EMConstantTerm;
  m_pfoCreatorSettings.m_hadConstantTerm = m_HadConstantTerm;
  
  // Calibration constants
  m_caloHitCreatorSettings.m_eCalToMip = m_ECalToMipCalibration; 
  m_caloHitCreatorSettings.m_hCalToMip = m_HCalToMipCalibration;
  m_caloHitCreatorSettings.m_eCalMipThreshold = m_ECalMipThreshold; 
  m_caloHitCreatorSettings.m_muonToMip = m_MuonToMipCalibration;
  m_caloHitCreatorSettings.m_hCalMipThreshold = m_HCalMipThreshold; 
  m_caloHitCreatorSettings.m_eCalToEMGeV = m_ECalToEMGeVCalibration;
  m_caloHitCreatorSettings.m_hCalToEMGeV = m_HCalToEMGeVCalibration;
  m_caloHitCreatorSettings.m_eCalToHadGeVEndCap = m_ECalToHadGeVCalibrationEndCap;
  m_caloHitCreatorSettings.m_eCalToHadGeVBarrel = m_ECalToHadGeVCalibrationBarrel; 
  m_caloHitCreatorSettings.m_hCalToHadGeV = m_HCalToHadGeVCalibration;
  m_caloHitCreatorSettings.m_muonDigitalHits = m_DigitalMuonHits; 
  m_caloHitCreatorSettings.m_muonHitEnergy = m_MuonHitEnergy; 
  m_caloHitCreatorSettings.m_maxHCalHitHadronicEnergy = m_MaxHCalHitHadronicEnergy; 
  m_caloHitCreatorSettings.m_nOuterSamplingLayers = m_NOuterSamplingLayers; 
  m_caloHitCreatorSettings.m_layersFromEdgeMaxRearDistance = m_LayersFromEdgeMaxRearDistance; 
  
  // Track relationship parameters
  m_trackCreatorSettings.m_shouldFormTrackRelationships = m_ShouldFormTrackRelationships; 
  // Initial track hit specifications
  m_trackCreatorSettings.m_minTrackHits = m_MinTrackHits;
  m_trackCreatorSettings.m_minFtdTrackHits = m_MinFtdTrackHits; 
  m_trackCreatorSettings.m_maxTrackHits = m_MaxTrackHits; 
  ////m_trackCreatorSettings.m_useOldTrackStateCalculation = m_UseOldTrackStateCalculation;
  // Track PFO usage parameters
  m_trackCreatorSettings.m_d0TrackCut = m_D0TrackCut; 
  m_trackCreatorSettings.m_z0TrackCut = m_Z0TrackCut;
  m_trackCreatorSettings.m_usingNonVertexTracks = m_UseNonVertexTracks;
  m_trackCreatorSettings.m_usingUnmatchedNonVertexTracks = m_UseUnmatchedNonVertexTracks;
  m_trackCreatorSettings.m_usingUnmatchedVertexTracks = m_UseUnmatchedVertexTracks;
  m_trackCreatorSettings.m_unmatchedVertexTrackMaxEnergy = m_UnmatchedVertexTrackMaxEnergy; 
  m_trackCreatorSettings.m_d0UnmatchedVertexTrackCut = m_D0UnmatchedVertexTrackCut; 
  m_trackCreatorSettings.m_z0UnmatchedVertexTrackCut = m_Z0UnmatchedVertexTrackCut;
  m_trackCreatorSettings.m_zCutForNonVertexTracks = m_ZCutForNonVertexTracks;
  // Track "reaches ecal" parameters
  m_trackCreatorSettings.m_reachesECalNTpcHits = m_ReachesECalNTpcHits;
  m_trackCreatorSettings.m_reachesECalNFtdHits = m_ReachesECalNFtdHits;
  m_trackCreatorSettings.m_reachesECalTpcOuterDistance = m_ReachesECalTpcOuterDistance;
  m_trackCreatorSettings.m_reachesECalMinFtdLayer = m_ReachesECalMinFtdLayer;
  m_trackCreatorSettings.m_reachesECalTpcZMaxDistance = m_ReachesECalTpcZMaxDistance;
  m_trackCreatorSettings.m_reachesECalFtdZMaxDistance = m_ReachesECalFtdZMaxDistance;
  m_trackCreatorSettings.m_curvatureToMomentumFactor = m_CurvatureToMomentumFactor;
  m_trackCreatorSettings.m_minTrackECalDistanceFromIp = m_MinTrackECalDistanceFromIp;
  // Final track quality parameters
  m_trackCreatorSettings.m_maxTrackSigmaPOverP = m_MaxTrackSigmaPOverP;
  m_trackCreatorSettings.m_minMomentumForTrackHitChecks = m_MinMomentumForTrackHitChecks;
  m_trackCreatorSettings.m_tpcMembraneMaxZ = m_TpcMembraneMaxZ;
  m_trackCreatorSettings.m_minTpcHitFractionOfExpected = m_MinTpcHitFractionOfExpected;
  m_trackCreatorSettings.m_minFtdHitsForTpcHitFraction = m_MinFtdHitsForTpcHitFraction;
  m_trackCreatorSettings.m_maxTpcInnerRDistance = m_MaxTpcInnerRDistance;
  
  
  // Additional geometry parameters
  m_geometryCreatorSettings.m_eCalEndCapInnerSymmetryOrder = m_ECalEndCapInnerSymmetryOrder;
  m_geometryCreatorSettings.m_eCalEndCapInnerPhiCoordinate = m_ECalEndCapInnerPhiCoordinate;
  m_geometryCreatorSettings.m_eCalEndCapOuterSymmetryOrder = m_ECalEndCapOuterSymmetryOrder;
  m_geometryCreatorSettings.m_eCalEndCapOuterPhiCoordinate = m_ECalEndCapOuterPhiCoordinate;
  m_geometryCreatorSettings.m_hCalEndCapInnerSymmetryOrder = m_HCalEndCapInnerSymmetryOrder;
  m_geometryCreatorSettings.m_hCalEndCapInnerPhiCoordinate = m_HCalEndCapInnerPhiCoordinate;
  m_geometryCreatorSettings.m_hCalEndCapOuterSymmetryOrder = m_HCalEndCapOuterSymmetryOrder;
  m_geometryCreatorSettings.m_hCalEndCapOuterPhiCoordinate = m_HCalEndCapOuterPhiCoordinate;
  m_geometryCreatorSettings.m_hCalRingInnerSymmetryOrder = m_HCalRingInnerSymmetryOrder;
  m_geometryCreatorSettings.m_hCalRingInnerPhiCoordinate = m_HCalRingInnerPhiCoordinate;
  m_geometryCreatorSettings.m_hCalRingOuterSymmetryOrder = m_HCalRingOuterSymmetryOrder; 
  m_geometryCreatorSettings.m_hCalRingOuterPhiCoordinate = m_HCalRingOuterPhiCoordinate;
  
  // For Strip Splitting method and also for hybrid ECAL
  m_caloHitCreatorSettings.m_stripSplittingOn = m_StripSplittingOn;
  m_caloHitCreatorSettings.m_useEcalScLayers = m_UseEcalScLayers;
  // Parameters for hybrid ECAL
  // Energy to MIP for Si-layers and Sc-layers, respectively.
  //Si
  m_caloHitCreatorSettings.m_eCalSiToMip = m_ECalSiToMipCalibration;
  //Sc
  m_caloHitCreatorSettings.m_eCalScToMip = m_ECalScToMipCalibration;
  // MipThreshold for Si-layers and Sc-layers, respectively.
  // Si
  m_caloHitCreatorSettings.m_eCalSiMipThreshold = m_ECalSiMipThreshold;
  //Sc
  m_caloHitCreatorSettings.m_eCalScMipThreshold = m_ECalScMipThreshold;
  // EcalToEM for Si-layers and Sc-layers, respectively.
  //Si
  m_caloHitCreatorSettings.m_eCalSiToEMGeV = m_ECalSiToEMGeVCalibration;
  //Sc
  m_caloHitCreatorSettings.m_eCalScToEMGeV = m_ECalScToEMGeVCalibration;
  // EcalToHad for Si-layers and Sc-layers of the endcaps, respectively.
  //Si
  m_caloHitCreatorSettings.m_eCalSiToHadGeVEndCap = m_ECalSiToHadGeVCalibrationEndCap;
  //Sc
  m_caloHitCreatorSettings.m_eCalScToHadGeVEndCap = m_ECalScToHadGeVCalibrationEndCap;
  // EcalToHad for Si-layers and Sc-layers of the barrel, respectively.
  //Si
  m_caloHitCreatorSettings.m_eCalSiToHadGeVBarrel = m_ECalSiToHadGeVCalibrationBarrel;
  //Sc
  m_caloHitCreatorSettings.m_eCalScToHadGeVBarrel = m_ECalScToHadGeVCalibrationBarrel;

  try
  {
      ISvcLocator* svcloc = serviceLocator();
      this->FinaliseSteeringParameters(svcloc);
      m_pPandora = new pandora::Pandora();
      m_pMCParticleCreator = new MCParticleCreator(m_mcParticleCreatorSettings, m_pPandora);
      m_pGeometryCreator = new GeometryCreator(m_geometryCreatorSettings, m_pPandora);
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pGeometryCreator->CreateGeometry(svcloc));
      m_pCaloHitCreator = new CaloHitCreator(m_caloHitCreatorSettings, m_pPandora, svcloc, 0);
      m_pTrackCreator = new TrackCreator(m_trackCreatorSettings, m_pPandora, svcloc);
      m_pPfoCreator = new PfoCreator(m_pfoCreatorSettings, m_pPandora);
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));
  }
  catch (pandora::StatusCodeException &statusCodeException)
  {
      std::cout << "Failed to initialize gaudi pandora: " << statusCodeException.ToString() << std::endl;
      throw statusCodeException;
  }
  catch (...)
  {
      std::cout << "Failed to initialize gaudi pandora: unrecognized exception" << std::endl;
      throw;
  }


  return GaudiAlgorithm::initialize();
}

StatusCode PandoraMatrixAlg::execute()
{
    
    try
    {
        std::cout<<"execute PandoraMatrixAlg"<<std::endl;
        
        updateMap();
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateMCParticles(*m_CollectionMaps));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pCaloHitCreator->CreateCaloHits(*m_CollectionMaps));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateCaloHitToMCParticleRelationships(*m_CollectionMaps, m_pCaloHitCreator->GetCalorimeterHitVector() ));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTrackAssociations(*m_CollectionMaps));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(*m_CollectionMaps));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pMCParticleCreator->CreateTrackToMCParticleRelationships(*m_CollectionMaps, m_pTrackCreator->GetTrackVector() ));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pPfoCreator->CreateParticleFlowObjects(*m_CollectionMaps, m_ClusterCollection_w, m_ReconstructedParticleCollection_w, m_VertexCollection_w));
        
        StatusCode sc0 = CreateMCRecoParticleAssociation();
        StatusCode sc = Ana();

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
        this->Reset();
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        std::cout << "Gaudi pandora failed to process event: " << statusCodeException.ToString() << std::endl;
        throw statusCodeException;
    }
    catch (...)
    {
        std::cout << "Gaudi pandora failed to process event: unrecognized exception" << std::endl;
        throw;
    }
  
  info() << "PandoraMatrixAlg Processed " << _nEvt << " events " << endmsg;
  _nEvt ++ ;

  return StatusCode::SUCCESS;
}

StatusCode PandoraMatrixAlg::finalize()
{
  info() << "Finalized. Processed " << _nEvt << " events " <<",saved tree with entries="<<m_tree->GetEntries()<< endmsg;
  m_fout->cd();
  m_tree->Write();
  m_fout->Close();
  delete m_pPandora;
  delete m_pGeometryCreator;
  delete m_pCaloHitCreator;
  delete m_pTrackCreator;
  delete m_pMCParticleCreator;
  delete m_pPfoCreator;
  return GaudiAlgorithm::finalize();
}




pandora::StatusCode PandoraMatrixAlg::RegisterUserComponents() const
{
    
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(*m_pPandora));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterBasicPlugins(*m_pPandora));

    
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterBFieldPlugin(*m_pPandora,
        m_settings.m_innerBField, m_settings.m_muonBarrelBField, m_settings.m_muonEndCapBField));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterNonLinearityEnergyCorrection(*m_pPandora,
        "NonLinearity", pandora::HADRONIC, m_settings.m_inputEnergyCorrectionPoints, m_settings.m_outputEnergyCorrectionPoints));
    
    

    return pandora::STATUS_CODE_SUCCESS;
}


void PandoraMatrixAlg::Reset()
{
    m_pCaloHitCreator->Reset();
    m_pTrackCreator->Reset();
    m_pMCParticleCreator->Reset();

    std::vector<int>()  .swap(m_pReco_PID   );
    std::vector< std::vector<int> >()  .swap(m_pReco_mc_id);
    std::vector< std::vector<float> >().swap(m_pReco_mc_weight);
    std::vector<float>().swap(m_pReco_mass);
    std::vector<float>().swap(m_pReco_energy);
    std::vector<float>().swap(m_pReco_px);
    std::vector<float>().swap(m_pReco_py);
    std::vector<float>().swap(m_pReco_pz);
    std::vector<float>().swap(m_pReco_charge);

    std::vector<int>()  .swap(m_mc_id);
    std::vector<int>()  .swap(m_mc_p_size);
    std::vector<int>()  .swap(m_mc_pid   );
    std::vector<float>().swap(m_mc_mass  );
    std::vector<float>().swap(m_mc_px    );
    std::vector<float>().swap(m_mc_py    );
    std::vector<float>().swap(m_mc_pz    );
    std::vector<float>().swap(m_mc_charge);

    std::vector<float>()  .swap(m_hits_x);
    std::vector<float>()  .swap(m_hits_y);
    std::vector<float>()  .swap(m_hits_z);
    std::vector<float>()  .swap(m_hits_E);

    m_hasConversion = 0;

    m_CollectionMaps->clear();
}

const pandora::Pandora *PandoraMatrixAlg::GetPandora() const
{
    if (NULL == m_pPandora)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_pPandora;
}
PandoraMatrixAlg::Settings::Settings() :
    m_innerBField(3.5f),
    m_muonBarrelBField(-1.5f),
    m_muonEndCapBField(0.01f)
{
}
CollectionMaps::CollectionMaps()
{
}
void CollectionMaps::clear()
{
collectionMap_MC.clear();
collectionMap_CaloHit.clear();
collectionMap_Vertex.clear();
collectionMap_Track.clear();
collectionMap_CaloRel.clear();
collectionMap_TrkRel.clear();
}



StatusCode PandoraMatrixAlg::updateMap()
{
    for(auto &v : m_dataHandles){
        try{
            if(m_collections[v.first]=="MCParticle"){
                auto handle = dynamic_cast<DataHandle<edm4hep::MCParticleCollection>*> (v.second);
                auto po = handle->get();
                if(po != NULL){
                    std::vector<edm4hep::MCParticle> v_mc;
                    m_CollectionMaps->collectionMap_MC [v.first] = v_mc;
                    for(unsigned int i=0 ; i< po->size(); i++) m_CollectionMaps->collectionMap_MC [v.first].push_back(po->at(i));
                    std::cout<<"saved col name="<<v.first<<std::endl;
                }
                else{
                std::cout<<"don't find col name="<<v.first<<std::endl;
                }
            }
            else if(m_collections[v.first]=="CalorimeterHit"){
                auto handle = dynamic_cast<DataHandle<edm4hep::CalorimeterHitCollection>*> (v.second);
                auto po = handle->get();
                if(po != NULL){
                    std::vector<edm4hep::CalorimeterHit> v_cal;
                    m_CollectionMaps->collectionMap_CaloHit[v.first] = v_cal ;
                    for(unsigned int i=0 ; i< po->size(); i++) m_CollectionMaps->collectionMap_CaloHit [v.first].push_back(po->at(i));
                    std::cout<<"saved col name="<<v.first<<std::endl;
                }
                else{
                std::cout<<"don't find col name="<<v.first<<std::endl;
                }
            }
            else if(m_collections[v.first]=="Track"){
                auto handle = dynamic_cast<DataHandle<edm4hep::TrackCollection>*> (v.second);
                auto po = handle->get();
                if(po != NULL){
                    std::vector<edm4hep::Track> v_cal;
                    m_CollectionMaps->collectionMap_Track[v.first] = v_cal ;
                    for(unsigned int i=0 ; i< po->size(); i++) m_CollectionMaps->collectionMap_Track [v.first].push_back(po->at(i));
                    std::cout<<"saved col name="<<v.first<<std::endl;
                }
                else{
                std::cout<<"don't find col name="<<v.first<<std::endl;
                }
            }
            else if(m_collections[v.first]=="Vertex"){
                auto handle = dynamic_cast<DataHandle<edm4hep::VertexCollection>*> (v.second);
                auto po = handle->get();
                if(po != NULL){
                    std::vector<edm4hep::Vertex> v_cal;
                    m_CollectionMaps->collectionMap_Vertex[v.first] = v_cal ;
                    for(unsigned int i=0 ; i< po->size(); i++) m_CollectionMaps->collectionMap_Vertex [v.first].push_back(po->at(i));
                    std::cout<<"saved col name="<<v.first<<std::endl;
                }
                else{
                std::cout<<"don't find col name="<<v.first<<std::endl;
                }
            }
            else if(m_collections[v.first]=="MCRecoCaloAssociation"){
                auto handle = dynamic_cast<DataHandle<edm4hep::MCRecoCaloAssociationCollection>*> (v.second);
                auto po = handle->get();
                if(po != NULL){
                    std::vector<edm4hep::MCRecoCaloAssociation> v_cal;
                    m_CollectionMaps->collectionMap_CaloRel[v.first] = v_cal ;
                    for(unsigned int i=0 ; i< po->size(); i++) m_CollectionMaps->collectionMap_CaloRel [v.first].push_back(po->at(i));
                    std::cout<<"saved col name="<<v.first<<std::endl;
                }
                else{
                std::cout<<"don't find col name="<<v.first<<std::endl;
                }
            }
            else if(m_collections[v.first]=="MCRecoTrackerAssociation"){
                auto handle = dynamic_cast<DataHandle<edm4hep::MCRecoTrackerAssociationCollection>*> (v.second);
                auto po = handle->get();
                if(po != NULL){
                    std::vector<edm4hep::MCRecoTrackerAssociation> v_cal;
                    m_CollectionMaps->collectionMap_TrkRel[v.first] = v_cal ;
                    for(unsigned int i=0 ; i< po->size(); i++) m_CollectionMaps->collectionMap_TrkRel [v.first].push_back(po->at(i));
                    std::cout<<"saved col name="<<v.first<<std::endl;
                }
                else{
                std::cout<<"don't find col name="<<v.first<<std::endl;
                }
            }
            else{
            std::cout<<"wrong type name for col :"<<v.first<<std::endl;
            }
        }//try
        catch(...){
            std::cout<<"don't find "<<v.first<<"in event"<<std::endl;
            std::cout<<"don't find  col name="<<v.first<<",with type="<<m_collections[v.first]<<" in this event"<<std::endl;
        }
    }
    return StatusCode::SUCCESS;
}


StatusCode PandoraMatrixAlg::Ana()
{
    int n_current = m_tree->GetEntries()+1;
    const edm4hep::ReconstructedParticleCollection* reco_col = m_ReconstructedParticleCollection_w.get();
    const edm4hep::MCRecoParticleAssociationCollection* reco_associa_col = m_MCRecoParticleAssociation_w.get();
    std::cout<<"reco_col size="<<reco_col->size()<<std::endl;
    for(int i=0; i<reco_col->size();i++)
    {
        std::cout<<"reco="<<i<<std::endl;
        const edm4hep::ReconstructedParticle pReco = reco_col->at(i);
        const float px = pReco.getMomentum()[0];
        const float py = pReco.getMomentum()[1];
        const float pz = pReco.getMomentum()[2];
        const float energy = pReco.getEnergy();
        const float mass = pReco.getMass();
        const float charge = pReco.getCharge();
        const int type = pReco.getType();
        std::cout<<"MYDBUG evt="<<n_current<<",rec i="<<i<<",particleId="<<type<<",mass="<<mass<<",charge="<<charge<<",energy="<<energy<<",px="<<px<<",py="<<py<<",pz="<<pz<<std::endl;
        m_pReco_PID.push_back(type);
        m_pReco_mass.push_back(mass);
        m_pReco_charge.push_back(charge);
        m_pReco_energy.push_back(energy);
        m_pReco_px.push_back(px);
        m_pReco_py.push_back(py);
        m_pReco_pz.push_back(pz);
        std::vector<int> tmp_mc_id;
        std::vector<float> tmp_mc_weight;
        for(int j=0; j < reco_associa_col->size(); j++)
        {
            if(reco_associa_col->at(j).getRec().id() != pReco.id() ) continue;
            std::cout<<"MC pid ="<<reco_associa_col->at(j).getSim().getPDG()<<",weight="<<reco_associa_col->at(j).getWeight()<<", px="<<reco_associa_col->at(j).getSim().getMomentum()[0]<<", py="<<reco_associa_col->at(j).getSim().getMomentum()[1]<<",pz="<<reco_associa_col->at(j).getSim().getMomentum()[2]<<std::endl;
            tmp_mc_id    .push_back(reco_associa_col->at(j).getSim().id());
            tmp_mc_weight.push_back(reco_associa_col->at(j).getWeight());
        }
        m_pReco_mc_id    .push_back(tmp_mc_id);
        m_pReco_mc_weight.push_back(tmp_mc_weight);
    }
    const edm4hep::MCParticleCollection*     MCParticle = nullptr;
    StatusCode sc = StatusCode::SUCCESS;
    sc =  getCol(m_mcParCol_r  , MCParticle );
    if (NULL != MCParticle   )  
    { 
        for(unsigned int i=0 ; i< MCParticle->size(); i++)
        {
            m_mc_id    .push_back(MCParticle->at(i).id());
            m_mc_p_size.push_back(MCParticle->at(i).parents_size());
            m_mc_pid   .push_back(MCParticle->at(i).getPDG());
            m_mc_mass  .push_back(MCParticle->at(i).getMass());
            m_mc_px    .push_back(MCParticle->at(i).getMomentum()[0]);
            m_mc_py    .push_back(MCParticle->at(i).getMomentum()[1]);
            m_mc_pz    .push_back(MCParticle->at(i).getMomentum()[2]);
            m_mc_charge.push_back(MCParticle->at(i).getCharge());
            //for(unsigned int j =0 ; j< MCParticle->at(i).daughters_size(); j++) da_pids.push_back( MCParticle->at(i).getDaughters(j).getPDG());
            if(MCParticle->at(i).parents_size()==0) std::cout<<"MYDBUG evt="<<n_current<<", mc i="<<i<<",px="<<MCParticle->at(i).getMomentum()[0]<<",py="<<MCParticle->at(i).getMomentum()[1]<<",pz="<<MCParticle->at(i).getMomentum()[2]<<std::endl;
            if (MCParticle->at(i).getPDG() != 22) continue;
            int hasEm = 0;
            int hasEp = 0;
            for(unsigned int j =0 ; j< MCParticle->at(i).daughters_size(); j++)
            {
                if      (MCParticle->at(i).getDaughters(j).getPDG() ==  11 ) hasEm=1;
                else if (MCParticle->at(i).getDaughters(j).getPDG() == -11 ) hasEp=1;
            }
            if(hasEm && hasEp) m_hasConversion=1;
        }
    }

    m_tree->Fill();
    return StatusCode::SUCCESS;
}

// create simple MCRecoParticleAssociation using calorimeter hit only
StatusCode PandoraMatrixAlg::CreateMCRecoParticleAssociation()
{
    edm4hep::MCRecoParticleAssociationCollection* pMCRecoParticleAssociationCollection  = m_MCRecoParticleAssociation_w.createAndPut();
    const edm4hep::ReconstructedParticleCollection* reco_col = m_ReconstructedParticleCollection_w.get();
    std::cout<<"CreateMCRecoParticleAssociation, reco_col size="<<reco_col->size()<<std::endl;
    for(int i=0; i<reco_col->size();i++)
    {
        std::map<int, edm4hep::ConstMCParticle> mc_map;
        std::map<int, float > id_edep_map;
        float tot_en = 0 ;
        const edm4hep::ReconstructedParticle pReco = reco_col->at(i);
        for(int j=0; j < pReco.clusters_size(); j++)
        {
            edm4hep::ConstCluster cluster = pReco.getClusters(j);
            for(int k=0; k < cluster.hits_size(); k++)
            {
                edm4hep::ConstCalorimeterHit hit = cluster.getHits(k);
                for(std::map<std::string, std::vector<edm4hep::MCRecoCaloAssociation> >::iterator iter = m_CollectionMaps->collectionMap_CaloRel.begin(); iter != m_CollectionMaps->collectionMap_CaloRel.end(); iter++)
                {
                    for(std::vector<edm4hep::MCRecoCaloAssociation>::iterator it = iter->second.begin(); it != iter->second.end(); it ++)
                    {
                        if(it->getRec().id() != hit.id()) continue;
                        for(std::vector<edm4hep::ConstCaloHitContribution>::const_iterator itc = it->getSim().contributions_begin(); itc != it->getSim().contributions_end(); itc++)
                        {
                            if(mc_map.find(itc->getParticle().id()) == mc_map.end()) mc_map[itc->getParticle().id()] = itc->getParticle() ;
                            if(id_edep_map.find(itc->getParticle().id()) != id_edep_map.end()) id_edep_map[itc->getParticle().id()] = id_edep_map[itc->getParticle().id()] + itc->getEnergy() ;
                            else                                                               id_edep_map[itc->getParticle().id()] = itc->getEnergy() ;
                            tot_en += itc->getEnergy() ;
                        }
                    }
                }
            }
        }
        for(std::map<int, edm4hep::ConstMCParticle>::iterator it = mc_map.begin(); it != mc_map.end(); it ++)
        {      
            edm4hep::MCRecoParticleAssociation association = pMCRecoParticleAssociationCollection->create();
            association.setRec(pReco);
            association.setSim(it->second);
            if(tot_en==0) 
            {
                association.setWeight(0);
                std::cout<<"Found 0 cluster energy"<<std::endl;
            }  
            else if(id_edep_map.find(it->first) != id_edep_map.end()) association.setWeight(id_edep_map[it->first]/tot_en);
            else std::cout<<"Error in creating MCRecoParticleAssociation"<<std::endl;
        }
    }
    return StatusCode::SUCCESS;
}

