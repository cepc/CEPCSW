#ifndef PandoraPFAlg_H
#define PandoraPFAlg_H

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCParticle.h" 
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCRecoCaloAssociation.h"
#include "edm4hep/MCRecoTrackerAssociation.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociation.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"

#include "Api/PandoraApi.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <iostream>
#include <random>
#include <string>
#include <unistd.h>


#include "CaloHitCreator.h"
#include "GeometryCreator.h"
#include "MCParticleCreator.h"
#include "PfoCreator.h"
#include "TrackCreator.h"


#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/MsgStream.h"

/* PandoraPFAlg ========== <br>
 * 
 */
namespace pandora {class Pandora;}


class CollectionMaps
{
public:
    CollectionMaps();
    void clear();
    std::map<std::string, std::vector<edm4hep::MCParticle> >     collectionMap_MC;
    std::map<std::string, std::vector<edm4hep::CalorimeterHit> > collectionMap_CaloHit;
    std::map<std::string, std::vector<edm4hep::Vertex> >         collectionMap_Vertex;
    std::map<std::string, std::vector<edm4hep::Track> >          collectionMap_Track;
    std::map<std::string, std::vector<edm4hep::MCRecoCaloAssociation> > collectionMap_CaloRel;
    std::map<std::string, std::vector<edm4hep::MCRecoTrackerAssociation> > collectionMap_TrkRel;
};



class PandoraPFAlg : public GaudiAlgorithm
{
 
public:
 
  PandoraPFAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;
 
  void FinaliseSteeringParameters(ISvcLocator* svcloc);
  pandora::StatusCode RegisterUserComponents() const;
  void Reset();
  typedef std::vector<float> FloatVector;
  typedef std::vector<std::string> StringVector;


  class Settings
  {
  public:
      /**
       *  @brief  Default constructor
       */
      Settings();

      std::string     m_pandoraSettingsXmlFile;           ///< The pandora settings xml file

      float           m_innerBField;                      ///< The bfield in the main tracker, ecal and hcal, units Tesla
      float           m_muonBarrelBField;                 ///< The bfield in the muon barrel, units Tesla
      float           m_muonEndCapBField;                 ///< The bfield in the muon endcap, units Tesla

      FloatVector     m_inputEnergyCorrectionPoints;      ///< The input energy points for non-linearity energy correction
      FloatVector     m_outputEnergyCorrectionPoints;     ///< The output energy points for non-linearity energy correction
  };


    /**
     *  @brief  Get address of the pandora instance
     * 
     *  @return address of the pandora instance
     */
    const pandora::Pandora *GetPandora() const;
    StatusCode updateMap();
    StatusCode Ana();
    StatusCode CreateMCRecoParticleAssociation();
protected:
 
  typedef std::vector<float> FloatVec;

  int _nEvt ;

 

  Gaudi::Property< std::string >              m_PandoraSettingsXmlFile { this, "PandoraSettingsDefault_xml", "PandoraSettingsDefault_wx.xml" };
  Gaudi::Property<int>                        m_NEventsToSkip                   { this, "NEventsToSkip", 0 };

  Gaudi::Property< std::vector<std::string> > m_TrackCollections{ this, "TrackCollections", {"MarlinTrkTracks"} };
  Gaudi::Property< std::vector<std::string> > m_ECalCaloHitCollections{ this, "ECalCaloHitCollections", {"ECALBarrel","ECALEndcap","ECALOther"} };
  Gaudi::Property< std::vector<std::string> > m_ECalReadOutNames      { this, "ECalReadOutNames", {"EcalBarrelCollection","EcalEndcapsCollection","ECALOther"} };
  Gaudi::Property< std::vector<std::string> > m_HCalCaloHitCollections{ this, "HCalCaloHitCollections", {"HCALBarrel","HCALEndcap","HCALOther"} };
  Gaudi::Property< std::vector<std::string> > m_HCalReadOutNames      { this, "HCalReadOutNames", {"HcalBarrelCollection","HcalEndcapsCollection","HCALOther"} };
  Gaudi::Property< std::vector<std::string> > m_LCalCaloHitCollections{ this, "LCalCaloHitCollections", {"LCAL"} };
  Gaudi::Property< std::vector<std::string> > m_LCalReadOutNames      { this, "LCalReadOutNames", {"LcalBarrelCollection"} };
  Gaudi::Property< std::vector<std::string> > m_LHCalCaloHitCollections{ this, "LHCalCaloHitCollections", {"LHCAL"} };
  Gaudi::Property< std::vector<std::string> > m_LHCalReadOutNames      { this, "LHCalReadOutNames", {"LHcalBarrelCollection"} };
  Gaudi::Property< std::vector<std::string> > m_MuonCaloHitCollections{ this, "MuonCaloHitCollections", {"MUON"} };
  Gaudi::Property< std::vector<std::string> > m_MuonCalReadOutNames      { this, "MuonCalReadOutNames", {"MuonCollection"} };
  Gaudi::Property< std::vector<std::string> > m_MCParticleCollections{ this, "MCParticleCollections", {"MCParticle"} };
  Gaudi::Property< std::vector<std::string> > m_RelCaloHitCollections{ this, "RelCaloHitCollections", {"RelationCaloHit","RelationMuonHit"} };
  Gaudi::Property< std::vector<std::string> > m_RelTrackCollections{ this, "RelTrackCollections", {"MarlinTrkTracksMCTruthLink"} };
  Gaudi::Property< std::vector<std::string> > m_KinkVertexCollections{ this, "KinkVertexCollections", {"KinkVertices"} };
  Gaudi::Property< std::vector<std::string> > m_ProngVertexCollections{ this, "ProngVertexCollections", {"ProngVertices"} };
  Gaudi::Property< std::vector<std::string> > m_SplitVertexCollections{ this, "SplitVertexCollections", {"SplitVertices"} };
  Gaudi::Property< std::vector<std::string> > m_V0VertexCollections{ this, "V0VertexCollections", {"V0Vertices"} };
  Gaudi::Property< std::string >              m_ClusterCollectionName { this, "ClusterCollectionName", "PandoraClusters" };
  Gaudi::Property< std::string >              m_PFOCollectionName { this, "PFOCollectionName", "PandoraPFOs" };
  Gaudi::Property<float>                      m_ECalToMipCalibration{ this, "ECalToMipCalibration", 160.0 };
  Gaudi::Property<float>                      m_HCalToMipCalibration{ this, "HCalToMipCalibration", 34.8 };
  Gaudi::Property<float>                      m_ECalMipThreshold{ this, "ECalMipThreshold", 0.5 };
  Gaudi::Property<float>                      m_HCalMipThreshold{ this, "HCalMipThreshold", 0.3 };
  Gaudi::Property<float>                      m_ECalToEMGeVCalibration{ this, "ECalToEMGeVCalibration", 1.007 };
  Gaudi::Property<float>                      m_HCalToEMGeVCalibration{ this, "HCalToEMGeVCalibration", 1.007 };
  Gaudi::Property<float>                      m_ECalToHadGeVCalibrationBarrel{ this, "ECalToHadGeVCalibrationBarrel", 1.12 };
  Gaudi::Property<float>                      m_ECalToHadGeVCalibrationEndCap{ this, "ECalToHadGeVCalibrationEndCap", 1.12 };
  Gaudi::Property<float>                      m_HCalToHadGeVCalibration{ this, "HCalToHadGeVCalibration", 1.07 };
  Gaudi::Property<float>                      m_MuonToMipCalibration{ this, "MuonToMipCalibration", 10.0 };
  Gaudi::Property<int>                        m_DigitalMuonHits{ this, "DigitalMuonHits", 0 };
  Gaudi::Property<float>                      m_MaxHCalHitHadronicEnergy{ this, "MaxHCalHitHadronicEnergy", 1.0 };
  Gaudi::Property<int>                        m_UseOldTrackStateCalculation{ this, "UseOldTrackStateCalculation", 0 };


  Gaudi::Property<float>                      m_AbsorberRadLengthECal{ this, "AbsorberRadLengthECal", 0.2854 };
  Gaudi::Property<float>                      m_AbsorberIntLengthECal{ this, "AbsorberIntLengthECal", 0.0101 };
  Gaudi::Property<float>                      m_AbsorberRadLengthHCal{ this, "AbsorberRadLengthHCal", 0.0569 };
  Gaudi::Property<float>                      m_AbsorberIntLengthHCal{ this, "AbsorberIntLengthHCal", 0.006  };
  Gaudi::Property<float>                      m_AbsorberRadLengthOther{ this, "AbsorberRadLengthOther", 0.0569  };
  Gaudi::Property<float>                      m_AbsorberIntLengthOther{ this, "AbsorberIntLengthOther", 0.006  };
  Gaudi::Property< std::string >              m_StartVertexCollectionName { this, "StartVertexCollectionName", "PandoraPFANewStartVertices" };
  Gaudi::Property< std::string >              m_StartVertexAlgorithmName { this, "StartVertexAlgorithmName", "PandoraPFANew" };
  Gaudi::Property<float>                      m_EMStochasticTerm{ this, "EMStochasticTerm", 0.17  };
  Gaudi::Property<float>                      m_HadStochasticTerm{ this, "HadStochasticTerm", 0.6  };
  Gaudi::Property<float>                      m_EMConstantTerm{ this, "EMConstantTerm", 0.01  };
  Gaudi::Property<float>                      m_HadConstantTerm{ this, "HadConstantTerm", 0.03  };
  Gaudi::Property<float>                      m_MuonHitEnergy{ this, "MuonHitEnergy", 0.5 };
  Gaudi::Property<int>                        m_NOuterSamplingLayers{ this, "NOuterSamplingLayers", 3 };
  Gaudi::Property<float>                      m_LayersFromEdgeMaxRearDistance{ this, "LayersFromEdgeMaxRearDistance", 250.f };
  Gaudi::Property<float>                      m_MuonBarrelBField{ this, "MuonBarrelBField", -1.5f };
  Gaudi::Property<float>                      m_MuonEndCapBField{ this, "MuonEndCapBField", 0.01f };
  Gaudi::Property<int>                        m_ShouldFormTrackRelationships{ this, "ShouldFormTrackRelationships", 1 };
  Gaudi::Property<int>                        m_MinTrackHits{ this, "MinTrackHits", 5 };
  Gaudi::Property<int>                        m_MinFtdTrackHits{ this, "MinFtdTrackHits", 0 };
  Gaudi::Property<int>                        m_MaxTrackHits{ this, "MaxTrackHits", 5000 };
  Gaudi::Property<float>                      m_D0TrackCut{ this, "D0TrackCut", 50. };
  Gaudi::Property<float>                      m_Z0TrackCut{ this, "Z0TrackCut", 50. };
  Gaudi::Property<int>                        m_UseNonVertexTracks{ this, "UseNonVertexTracks", 1 };
  Gaudi::Property<int>                        m_UseUnmatchedNonVertexTracks{ this, "UseUnmatchedNonVertexTracks", 0 };
  Gaudi::Property<int>                        m_UseUnmatchedVertexTracks{ this, "UseUnmatchedVertexTracks", 1 };
  Gaudi::Property<float>                      m_UnmatchedVertexTrackMaxEnergy{ this, "UnmatchedVertexTrackMaxEnergy", 5. };
  Gaudi::Property<float>                      m_D0UnmatchedVertexTrackCut{ this, "D0UnmatchedVertexTrackCut", 5. };
  Gaudi::Property<float>                      m_Z0UnmatchedVertexTrackCut{ this, "Z0UnmatchedVertexTrackCut", 5. };
  Gaudi::Property<float>                      m_ZCutForNonVertexTracks{ this, "ZCutForNonVertexTracks", 250. };
  Gaudi::Property<int>                        m_ReachesECalNTpcHits{ this, "ReachesECalNTpcHits", 11 };
  Gaudi::Property<int>                        m_ReachesECalNFtdHits{ this, "ReachesECalNFtdHits", 4 };
  Gaudi::Property<float>                      m_ReachesECalTpcOuterDistance{ this, "ReachesECalTpcOuterDistance", -100. };
  Gaudi::Property<int>                        m_ReachesECalMinFtdLayer{ this, "ReachesECalMinFtdLayer", 9 };
  Gaudi::Property<float>                      m_ReachesECalTpcZMaxDistance{ this, "ReachesECalTpcZMaxDistance", -50. };
  Gaudi::Property<float>                      m_ReachesECalFtdZMaxDistance{ this, "ReachesECalFtdZMaxDistance", -1. };
  Gaudi::Property<float>                      m_CurvatureToMomentumFactor{ this, "CurvatureToMomentumFactor", 0.3 / 2000. };
  Gaudi::Property<float>                      m_MinTrackECalDistanceFromIp{ this, "MinTrackECalDistanceFromIp", 100. };

  Gaudi::Property<float>                      m_MaxTrackSigmaPOverP             { this, "MaxTrackSigmaPOverP", 0.15 };
  Gaudi::Property<float>                      m_MinMomentumForTrackHitChecks    { this, "MinMomentumForTrackHitChecks", 1. };
  Gaudi::Property<float>                      m_TpcMembraneMaxZ                 { this, "TpcMembraneMaxZ", 10. };
  Gaudi::Property<float>                      m_MinTpcHitFractionOfExpected     { this, "MinTpcHitFractionOfExpected", 0.20 };
  Gaudi::Property<int>                        m_MinFtdHitsForTpcHitFraction     { this, "MinFtdHitsForTpcHitFraction", 2 };
  Gaudi::Property<float>                      m_MaxTpcInnerRDistance            { this, "MaxTpcInnerRDistance", 50.0 };
  Gaudi::Property<int>                        m_ECalEndCapInnerSymmetryOrder    { this, "ECalEndCapInnerSymmetryOrder", 4 };
  Gaudi::Property<float>                      m_ECalEndCapInnerPhiCoordinate    { this, "ECalEndCapInnerPhiCoordinate", 0. };
  Gaudi::Property<int>                        m_ECalEndCapOuterSymmetryOrder    { this, "ECalEndCapOuterSymmetryOrder", 8 };
  Gaudi::Property<float>                      m_ECalEndCapOuterPhiCoordinate    { this, "ECalEndCapOuterPhiCoordinate", 0. };
  Gaudi::Property<int>                        m_HCalEndCapInnerSymmetryOrder    { this, "HCalEndCapInnerSymmetryOrder", 4 };
  Gaudi::Property<float>                      m_HCalEndCapInnerPhiCoordinate    { this, "HCalEndCapInnerPhiCoordinate", 0. };
  Gaudi::Property<int>                        m_HCalEndCapOuterSymmetryOrder    { this, "HCalEndCapOuterSymmetryOrder", 16 };
  Gaudi::Property<float>                      m_HCalEndCapOuterPhiCoordinate    { this, "HCalEndCapOuterPhiCoordinate", 0. };
  Gaudi::Property<int>                        m_HCalRingInnerSymmetryOrder      { this, "HCalRingInnerSymmetryOrder",   8  };
  Gaudi::Property<float>                      m_HCalRingInnerPhiCoordinate      { this, "HCalRingInnerPhiCoordinate",   0. };
  Gaudi::Property<int>                        m_HCalRingOuterSymmetryOrder      { this, "HCalRingOuterSymmetryOrder",   16 };
  Gaudi::Property<float>                      m_HCalRingOuterPhiCoordinate      { this, "HCalRingOuterPhiCoordinate",   0. };
  Gaudi::Property<bool>                       m_StripSplittingOn                { this, "StripSplittingOn", false };
  Gaudi::Property<bool>                       m_UseEcalScLayers                 { this, "UseEcalScLayers", false };
  Gaudi::Property<float>                      m_ECalSiToMipCalibration          { this, "ECalSiToMipCalibration", 1. };
  Gaudi::Property<float>                      m_ECalScToMipCalibration          { this, "ECalScToMipCalibration", 1. };
  Gaudi::Property<float>                      m_ECalSiMipThreshold              { this, "ECalSiMipThreshold", 0. };
  Gaudi::Property<float>                      m_ECalScMipThreshold              { this, "ECalScMipThreshold", 0. };
  Gaudi::Property<float>                      m_ECalSiToEMGeVCalibration        { this, "ECalSiToEMGeVCalibration", 1. };
  Gaudi::Property<float>                      m_ECalScToEMGeVCalibration        { this, "ECalScToEMGeVCalibration", 1. };
  Gaudi::Property<float>                      m_ECalSiToHadGeVCalibrationEndCap { this, "ECalSiToHadGeVCalibrationEndCap", 1. };
  Gaudi::Property<float>                      m_ECalScToHadGeVCalibrationEndCap { this, "ECalScToHadGeVCalibrationEndCap", 1. };
  Gaudi::Property<float>                      m_ECalSiToHadGeVCalibrationBarrel { this, "ECalSiToHadGeVCalibrationBarrel", 1. };
  Gaudi::Property<float>                      m_ECalScToHadGeVCalibrationBarrel { this, "ECalScToHadGeVCalibrationBarrel", 1. };

  Gaudi::Property<FloatVector>                m_InputEnergyCorrectionPoints { this, "InputEnergyCorrectionPoints", {} };
  Gaudi::Property<FloatVector>                m_OutputEnergyCorrectionPoints { this, "OutputEnergyCorrectionPoints", {} };


  static pandora::Pandora        *m_pPandora;
  GeometryCreator                *m_pGeometryCreator;             ///< The geometry creator
  CaloHitCreator                 *m_pCaloHitCreator;              ///< The calo hit creator
  TrackCreator                   *m_pTrackCreator;                ///< The track creator
  MCParticleCreator              *m_pMCParticleCreator;           ///< The mc particle creator
  PfoCreator                     *m_pPfoCreator;                  ///< The pfo creator
 
  Settings                        m_settings;                     ///< The settings for the pandora pfa new algo
  CollectionMaps                  *m_CollectionMaps;               ///< The settings for the pandora pfa new algo
  GeometryCreator::Settings       m_geometryCreatorSettings;      ///< The geometry creator settings
  TrackCreator::Settings          m_trackCreatorSettings;         ///< The track creator settings
  CaloHitCreator::Settings        m_caloHitCreatorSettings;       ///< The calo hit creator settings
  MCParticleCreator::Settings     m_mcParticleCreatorSettings;    ///< The mc particle creator settings
  PfoCreator::Settings            m_pfoCreatorSettings;           ///< The pfo creator settings

  std::string                     m_detectorName;                 ///< The detector name
  unsigned int                    m_nRun;                         ///< The run number
  unsigned int                    m_nEvent;                       ///< The event number
  //### For Ana #################
  NTuple::Tuple* m_tuple = nullptr ;
  NTuple::Item<long>   m_n_mc;
  NTuple::Item<long>   m_n_rec;
  NTuple::Item<int>   m_hasConversion;
  NTuple::Item<int>   m_marlinTrack;
  NTuple::Array<int  > m_pReco_PID;    
  NTuple::Array<float> m_pReco_mass;
  NTuple::Array<float> m_pReco_energy;
  NTuple::Array<float> m_pReco_px;
  NTuple::Array<float> m_pReco_py;
  NTuple::Array<float> m_pReco_pz;
  NTuple::Array<float> m_pReco_charge;
  NTuple::Array<int>   m_mc_p_size;
  NTuple::Array<int>   m_mc_pid   ;
  NTuple::Array<float> m_mc_mass  ;
  NTuple::Array<float> m_mc_px    ;
  NTuple::Array<float> m_mc_py    ;
  NTuple::Array<float> m_mc_pz    ;
  NTuple::Array<float> m_mc_charge;


  Gaudi::Property<int> m_max_mc {this, "max_mc", 1000,""};
  Gaudi::Property<int> m_max_rec {this, "max_rec", 1000,""};

  
  Gaudi::Property<bool> m_debug {this, "debug", false,"if do debug"};
  Gaudi::Property<bool> m_WriteAna {this, "WriteAna", false,"if do ana"};
  Gaudi::Property<bool> m_use_dd4hep_geo{this, "use_dd4hep_geo", false,"choose if use geo info from dd4hep"};
  Gaudi::Property<bool> m_use_dd4hep_decoder {this, "use_dd4hep_decoder", true,"if use decoder from dd4hep"};
  Gaudi::Property<bool> m_use_preshower {this, "use_preshower", false,"if use preshower layer for calorimeter"};
  //######################
  std::map< std::string, std::string > m_collections;
  Gaudi::Property<std::vector<std::string>> m_readCols{this, "collections", {}, "Places of collections to read"};
 //the map of collection name to its corresponding DataHandle
  std::map<std::string, DataObjectHandleBase*> m_dataHandles;
  
  DataHandle<edm4hep::MCParticleCollection>     m_mcParCol_r  {"MCParticle", Gaudi::DataHandle::Reader, this};

  DataHandle<edm4hep::ClusterCollection>                m_ClusterCollection_w {"PandoraClusters",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::ReconstructedParticleCollection>  m_ReconstructedParticleCollection_w {"PandoraPFOs"    ,Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::VertexCollection>                 m_VertexCollection_w {"PandoraPFANewStartVertices",Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoParticleAssociationCollection>  m_MCRecoParticleAssociation_w {"pfoMCRecoParticleAssociation",Gaudi::DataHandle::Writer, this};

};

#endif
