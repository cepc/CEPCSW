/**
 * 
 *  @brief  Header file for the calo hit creator class.
 * 
 *  $Log: $
 */

#ifndef CALO_HIT_CREATOR_H
#define CALO_HIT_CREATOR_H 1

#include "GaudiKernel/ISvcLocator.h"
#include "edm4hep/CalorimeterHit.h"

#include "gear/LayerLayout.h"

#include "Api/PandoraApi.h"

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DD4hep/BitFieldCoder.h"


#include <string>

typedef std::vector<edm4hep::CalorimeterHit *> CalorimeterHitVector;

namespace gear { class GearMgr; }

class CollectionMaps;
/**
 *  @brief  CaloHitCreator class
 */
class CaloHitCreator
{
public:
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        StringVector    m_eCalCaloHitCollections;               ///< The ecal calorimeter hit collections
        StringVector    m_hCalCaloHitCollections;               ///< The hcal calorimeter hit collections
        StringVector    m_lCalCaloHitCollections;               ///< The lcal calorimeter hit collections
        StringVector    m_lHCalCaloHitCollections;              ///< The lhcal calorimeter hit collections
        StringVector    m_muonCaloHitCollections;               ///< The muon calorimeter hit collections

        float           m_absorberRadLengthECal;                ///< The absorber radiation length in the ECal
        float           m_absorberIntLengthECal;                ///< The absorber interaction length in the ECal
        float           m_absorberRadLengthHCal;                ///< The absorber radiation length in the HCal
        float           m_absorberIntLengthHCal;                ///< The absorber interaction length in the HCal
        float           m_absorberRadLengthOther;               ///< The absorber radiation length in other detector regions
        float           m_absorberIntLengthOther;               ///< The absorber interaction length in other detector regions

        float           m_eCalToMip;                            ///< The calibration from deposited ECal energy to mip
        float           m_hCalToMip;                            ///< The calibration from deposited HCal energy to mip
        float           m_muonToMip;                            ///< The calibration from deposited Muon energy to mip
        float           m_eCalMipThreshold;                     ///< Threshold for creating calo hits in the ECal, units mip
        float           m_hCalMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip
        float           m_muonMipThreshold;                     ///< Threshold for creating calo hits in the HCal, units mip

        float           m_eCalToEMGeV;                          ///< The calibration from deposited ECal energy to EM energy
        float           m_eCalToHadGeVBarrel;                   ///< The calibration from deposited ECal barrel energy to hadronic energy
        float           m_eCalToHadGeVEndCap;                   ///< The calibration from deposited ECal endcap energy to hadronic energy
        float           m_hCalToEMGeV;                          ///< The calibration from deposited HCal energy to EM energy
        float           m_hCalToHadGeV;                         ///< The calibration from deposited HCal energy to hadronic energy
        int             m_muonDigitalHits;                      ///< Muon hits are treated as digital (energy from hit count)
        float           m_muonHitEnergy;                        ///< The energy for a digital muon calorimeter hit, units GeV

        float           m_maxHCalHitHadronicEnergy;             ///< The maximum hadronic energy allowed for a single hcal hit
        int             m_nOuterSamplingLayers;                 ///< Number of layers from edge for hit to be flagged as an outer layer hit
        float           m_layersFromEdgeMaxRearDistance;        ///< Maximum number of layers from candidate outer layer hit to rear of detector

        int             m_hCalEndCapInnerSymmetryOrder;         ///< HCal end cap inner symmetry order (missing from ILD00 gear file)
        float           m_hCalEndCapInnerPhiCoordinate;         ///< HCal end cap inner phi coordinate (missing from ILD00 gear file)

        // For Strip Splitting method and hybrid ECAL.
        int             m_stripSplittingOn;                     ///< To use SSA, this should be true (default is false)
        int             m_useEcalScLayers;                      ///< To use scintillator layers ~ hybrid ECAL, this should be true (default is false)
        float           m_eCalSiToMip;                          ///< The calibration from deposited Si-layer energy to mip
        float           m_eCalScToMip;                          ///< The calibration from deposited Sc-layer energy to mip
        float           m_eCalSiMipThreshold;                   ///< Threshold for creating calo hits in the Si-layers of ECAL, units mip
        float           m_eCalScMipThreshold;                   ///< Threshold for creating calo hits in the Sc-layers of ECAL, units mip
        float           m_eCalSiToEMGeV;                        ///< The calibration from deposited Si-layer energy to EM energy
        float           m_eCalScToEMGeV;                        ///< The calibration from deposited Sc-layer energy to EM energy
        float           m_eCalSiToHadGeVBarrel;                 ///< The calibration from deposited Si-layer energy on the enecaps to hadronic energy
        float           m_eCalScToHadGeVBarrel;                 ///< The calibration from deposited Sc-layer energy on the endcaps to hadronic energy
        float           m_eCalSiToHadGeVEndCap;                 ///< The calibration from deposited Si-layer energy on the enecaps to hadronic energy
        float           m_eCalScToHadGeVEndCap;                 ///< The calibration from deposited Sc-layer energy on the endcaps to hadronic energy
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora, ISvcLocator* svcloc, bool encoder_style);

    /**
     *  @brief  Destructor
     */
     ~CaloHitCreator();

    /**
     *  @brief  Create calo hits
     * 
     */    
    pandora::StatusCode CreateCaloHits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Get the calorimeter hit vector
     * 
     *  @return The calorimeter hit vector
     */
    const CalorimeterHitVector &GetCalorimeterHitVector() const;

    /**
     *  @brief  Reset the calo hit creator
     */
    void Reset();

private:
    /**
     *  @brief  Create ecal calo hits
     * 
     */
    pandora::StatusCode CreateECalCaloHits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Create hcal calo hits
     * 
     */
    pandora::StatusCode CreateHCalCaloHits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Create muon calo hits
     * 
     */
    pandora::StatusCode CreateMuonCaloHits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Create lcal calo hits
     * 
     */    
    pandora::StatusCode CreateLCalCaloHits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Create lhcal calo hits
     * 
     */
    pandora::StatusCode CreateLHCalCaloHits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Get common calo hit properties: position, parent address, input energy and time
     * 
     */
    void GetCommonCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

    /**
     *  @brief  Get end cap specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layerLayout the gear end cap layer layout
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
    void GetEndCapCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
        PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const;

    /**
     *  @brief  Get barrel specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  layerLayout the gear barrel layer layout
     *  @param  barrelSymmetryOrder the barrel order of symmetry
     *  @param  barrelPhi0 the barrel orientation
     *  @param  staveNumber the stave number
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
    void GetBarrelCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
        unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
        float &absorberCorrection) const;

    void GetBarrelCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer> &layerLayout,
        unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
        float &absorberCorrection) const;
    /**
     *  @brief  Get number of active layers from position of a calo hit to the edge of the detector
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     */
    int GetNLayersFromEdge(const edm4hep::CalorimeterHit *const pCaloHit) const;

    /**
     *  @brief  Get the maximum radius of a calo hit in a polygonal detector structure
     * 
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  symmetryOrder the symmetry order
     *  @param  phi0 the angular orientation
     * 
     *  @return the maximum radius
     */
    float GetMaximumRadius(const edm4hep::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const;

    void GetCoding(const edm4hep::CalorimeterHit* pCaloHit, long& sys, long& x, long& y, long& z) const ;
    int GetBarrelLayer(const edm4hep::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer> &layerLayout) const;
    /**
     *  @brief  Get the layer coding string from the provided cell id encoding string
     * 
     *  @param  encodingString the cell id encoding string
     * 
     *  @return the layer coding string
     */
    std::string GetLayerCoding(const std::string &encodingString) const;

    /**
     *  @brief  Get the stave coding string from the provided cell id encoding string
     * 
     *  @param  encodingString the cell id encoding string
     * 
     *  @return the stave coding string
     */
    std::string GetStaveCoding(const std::string &encodingString) const;

    const Settings                      m_settings;                         ///< The calo hit creator settings

    const pandora::Pandora             *m_pPandora;                         ///< Address of the pandora object to create calo hits

    float                               m_eCalBarrelOuterZ;                 ///< ECal barrel outer z coordinate
    float                               m_hCalBarrelOuterZ;                 ///< HCal barrel outer z coordinate
    float                               m_muonBarrelOuterZ;                 ///< Muon barrel outer z coordinate
    float                               m_coilOuterR;                       ///< Coil outer r coordinate

    float                               m_eCalBarrelInnerPhi0;              ///< ECal barrel inner phi0 coordinate
    unsigned int                        m_eCalBarrelInnerSymmetry;          ///< ECal barrel inner symmetry order
    float                               m_hCalBarrelInnerPhi0;              ///< HCal barrel inner phi0 coordinate
    unsigned int                        m_hCalBarrelInnerSymmetry;          ///< HCal barrel inner symmetry order
    float                               m_muonBarrelInnerPhi0;              ///< Muon barrel inner phi0 coordinate
    unsigned int                        m_muonBarrelInnerSymmetry;          ///< Muon barrel inner symmetry order

    float                               m_hCalEndCapOuterR;                 ///< HCal endcap outer r coordinate
    float                               m_hCalEndCapOuterZ;                 ///< HCal endcap outer z coordinate
    float                               m_hCalBarrelOuterR;                 ///< HCal barrel outer r coordinate
    float                               m_hCalBarrelOuterPhi0;              ///< HCal barrel outer phi0 coordinate
    unsigned int                        m_hCalBarrelOuterSymmetry;          ///< HCal barrel outer symmetry order

    float                               m_hCalBarrelLayerThickness;         ///< HCal barrel layer thickness
    float                               m_hCalEndCapLayerThickness;         ///< HCal endcap layer thickness

    CalorimeterHitVector                m_calorimeterHitVector;             ///< The calorimeter hit vector
    std::string                         m_encoder_str;
    std::string                         m_encoder_str_MUON ; 
    std::string                         m_encoder_str_LCal ; 
    std::string                         m_encoder_str_LHCal; 
    gear::GearMgr* _GEAR;
    dd4hep::Detector* m_dd4hep;
    dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
    std::string                         m_compact ;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CalorimeterHitVector &CaloHitCreator::GetCalorimeterHitVector() const
{
    return m_calorimeterHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void CaloHitCreator::Reset()
{
    m_calorimeterHitVector.clear();
}

#endif // #ifndef CALO_HIT_CREATOR_H
