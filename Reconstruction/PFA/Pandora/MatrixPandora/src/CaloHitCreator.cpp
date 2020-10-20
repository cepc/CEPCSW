/**
 * 
 *  @brief  Implementation of the calo hit creator class.
 * 
 *  $Log: $
 */


#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/GearDistanceProperties.h"
#include "gear/GearPointProperties.h"
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"

#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Objects.h>

#include "UTIL/CellIDDecoder.h"
#include "cellIDDecoder.h"
#include "GaudiKernel/IService.h"
#include "GearSvc/IGearSvc.h"

#include "PandoraMatrixAlg.h"
#include "CaloHitCreator.h"

#include <algorithm>
#include <cmath>
#include <limits>


CaloHitCreator::CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora, ISvcLocator* svcloc, bool encoder_style) :
    m_settings(settings),
    m_pPandora(pPandora)
{
    m_encoder_str = ""; 
    m_encoder_str_MUON = ""; 
    m_encoder_str_LCal = ""; 
    m_encoder_str_LHCal = ""; 
    if(encoder_style==0) 
    {
        m_encoder_str     = "M:3,S-1:3,I:9,J:9,K-1:6";// LCIO style
        m_encoder_str_MUON="S-1:4,M:3,K-1:6,I:16,GRZone:3,J:32:16";
        m_encoder_str_LCal="I:10,J:10,K:10,S-1:2";
        m_encoder_str_LHCal=m_encoder_str;
    }
    IGearSvc*  iSvc = 0;
    StatusCode sc = svcloc->service("GearSvc", iSvc, false);
    if ( !sc ) 
    {
        throw "Failed to find GearSvc ...";
    }
    _GEAR = iSvc->getGearMgr();


    IGeomSvc*  Svc = 0;
    StatusCode sc1 = svcloc->service("GeomSvc", Svc, false);
    if ( !sc1 )
    {
        throw "Failed to find GeomSvc ...";
    }
    m_dd4hep = Svc->lcdd();
    const dd4hep::DetElement &detElement = m_dd4hep->detector("CaloDetector");
    dd4hep::rec::LayeredCalorimeterData* Data = detElement.extension<dd4hep::rec::LayeredCalorimeterData>() ;
    if(!Data) throw "Failed to get LayeredCalorimeterData ...";
    m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);

    m_eCalBarrelOuterZ        = Data->extent[3];
    m_eCalBarrelInnerPhi0     = Data->inner_phi0;
    m_eCalBarrelInnerSymmetry = Data->inner_symmetry;
    m_hCalBarrelOuterZ        = (_GEAR->getHcalBarrelParameters().getExtent()[3]);
    m_muonBarrelOuterZ        = (_GEAR->getYokeBarrelParameters().getExtent()[3]);
    m_coilOuterR              = (_GEAR->getGearParameters("CoilParameters").getDoubleVal("Coil_cryostat_outer_radius"));
    m_hCalBarrelInnerPhi0     = (_GEAR->getHcalBarrelParameters().getPhi0());
    m_hCalBarrelInnerSymmetry = (_GEAR->getHcalBarrelParameters().getSymmetryOrder());
    m_muonBarrelInnerPhi0     = (_GEAR->getYokeBarrelParameters().getPhi0());
    m_muonBarrelInnerSymmetry = (_GEAR->getYokeBarrelParameters().getSymmetryOrder());
    m_hCalEndCapOuterR        = (_GEAR->getHcalEndcapParameters().getExtent()[1]);
    m_hCalEndCapOuterZ        = (_GEAR->getHcalEndcapParameters().getExtent()[3]);
    m_hCalBarrelOuterR        = (_GEAR->getHcalBarrelParameters().getExtent()[1]);
    m_hCalBarrelOuterPhi0     =((std::find(_GEAR->getHcalBarrelParameters().getIntKeys().begin(),
        _GEAR->getHcalBarrelParameters().getIntKeys().end(),
        "Hcal_outer_polygon_phi0") != _GEAR->getHcalBarrelParameters().getIntKeys().end() ?
        _GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0")
        : 0));
    m_hCalBarrelOuterSymmetry = ((std::find(_GEAR->getHcalBarrelParameters().getIntKeys().begin(),
        _GEAR->getHcalBarrelParameters().getIntKeys().end(),
        "Hcal_outer_polygon_order") != _GEAR->getHcalBarrelParameters().getIntKeys().end() ?
        _GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order")
        : 0));


    const gear::LayerLayout &hCalEndCapLayerLayout(_GEAR->getHcalEndcapParameters().getLayerLayout());
    const gear::LayerLayout &hCalBarrelLayerLayout(_GEAR->getHcalBarrelParameters().getLayerLayout()); 
    m_hCalEndCapLayerThickness = hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1);
    m_hCalBarrelLayerThickness = hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1);
    if ((m_hCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()) || (m_hCalBarrelLayerThickness < std::numeric_limits<float>::epsilon()))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);


}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitCreator::~CaloHitCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateCaloHits(const CollectionMaps& collectionMaps)
{
    
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits (collectionMaps));
    // not used for ECAL Matrix
    /*
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits (collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits (collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLCalCaloHits (collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLHCalCaloHits(collectionMaps));
    */
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateECalCaloHits(const CollectionMaps& collectionMaps)
{
    std::cout<<"start CreateECalCaloHits:"<<std::endl;
    for (StringVector::const_iterator iter = m_settings.m_eCalCaloHitCollections.begin(), iterEnd = m_settings.m_eCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

            std::cout<<(*iter)<<"has nElements="<<nElements<<std::endl;
            if (0 == nElements)
                continue;

            const gear::LayerLayout &endcapLayerLayout(_GEAR->getEcalEndcapParameters().getLayerLayout());
            const dd4hep::DetElement &detElement = m_dd4hep->detector("CaloDetector");
            dd4hep::rec::LayeredCalorimeterData* Data = detElement.extension<dd4hep::rec::LayeredCalorimeterData>() ;
            const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>& barrelLayerLayout = Data->layers;
            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str);
            const std::string layerCodingString(m_encoder_str);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));
            const std::string staveCoding(this->GetStaveCoding(layerCodingString));

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    const edm4hep::CalorimeterHit& pCaloHit0 = pCaloHitCollection.at(i);
                    const edm4hep::CalorimeterHit* pCaloHit = &(pCaloHit0);

                    if (NULL == pCaloHit)
                        throw ("CreateECalCaloHits pCaloHit Collection type mismatch");

                    float eCalToMip(m_settings.m_eCalToMip), eCalMipThreshold(m_settings.m_eCalMipThreshold), eCalToEMGeV(m_settings.m_eCalToEMGeV),
                        eCalToHadGeVBarrel(m_settings.m_eCalToHadGeVBarrel), eCalToHadGeVEndCap(m_settings.m_eCalToHadGeVEndCap);

                    // Hybrid ECAL including pure ScECAL.
                    if (m_settings.m_useEcalScLayers)
                    {
                        std::string collectionName(*iter);
                        std::transform(collectionName.begin(), collectionName.end(), collectionName.begin(), ::tolower);

                        if (collectionName.find("ecal", 0) == std::string::npos)
                            std::cout << "WARNING: mismatching hybrid Ecal collection name. " << collectionName << std::endl;

                        if (collectionName.find("si", 0) != std::string::npos)
                        {
                             eCalToMip = m_settings.m_eCalSiToMip;
                             eCalMipThreshold = m_settings.m_eCalSiMipThreshold;
                             eCalToEMGeV = m_settings.m_eCalSiToEMGeV;
                             eCalToHadGeVBarrel = m_settings.m_eCalSiToHadGeVBarrel;
                             eCalToHadGeVEndCap = m_settings.m_eCalSiToHadGeVEndCap;
                        }
                        else if (collectionName.find("sc", 0) != std::string::npos)
                        {
                             eCalToMip = m_settings.m_eCalScToMip;
                             eCalMipThreshold = m_settings.m_eCalScMipThreshold;
                             eCalToEMGeV = m_settings.m_eCalScToEMGeV;
                             eCalToHadGeVBarrel = m_settings.m_eCalScToHadGeVBarrel;
                             eCalToHadGeVEndCap = m_settings.m_eCalScToHadGeVEndCap;
                        }
                    }

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_eCalBarrelOuterZ)
                    {
                        caloHitParameters.m_layer = GetBarrelLayer(pCaloHit, barrelLayerLayout);
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_eCalBarrelInnerSymmetry, m_eCalBarrelInnerPhi0, 6, caloHitParameters, absorberCorrection);//6
                        caloHitParameters.m_hadronicEnergy = eCalToHadGeVBarrel * pCaloHit->getEnergy();
                        
                    }
                    else
                    { // will not be used for ECAL Matrix
                        caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()] + 1;
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                        caloHitParameters.m_hadronicEnergy = eCalToHadGeVEndCap * pCaloHit->getEnergy();
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = eCalToEMGeV * pCaloHit->getEnergy();

                    // ATTN If using strip splitting, must correct cell sizes for use in PFA to minimum of strip width and strip length
                    if (m_settings.m_stripSplittingOn)
                    {
                        const float splitCellSize(std::min(caloHitParameters.m_cellSize0.Get(), caloHitParameters.m_cellSize1.Get()));
                        caloHitParameters.m_cellSize0 = splitCellSize;
                        caloHitParameters.m_cellSize1 = splitCellSize;
                    }

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(const_cast<edm4hep::CalorimeterHit*>(pCaloHit));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout<<"ERROR Failed to extract ecal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout<<"WARNING Failed to extract ecal calo hit: " <<  std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout<< "Failed to extract ecal calo hit collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateHCalCaloHits(const CollectionMaps& collectionMaps)
{
    for (StringVector::const_iterator iter = m_settings.m_hCalCaloHitCollections.begin(), iterEnd = m_settings.m_hCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

            std::cout<<(*iter)<<"has nElements="<<nElements<<std::endl;
            if (0 == nElements)
                continue;

            const gear::LayerLayout &endcapLayerLayout(_GEAR->getHcalEndcapParameters().getLayerLayout());
            const gear::LayerLayout &barrelLayerLayout(_GEAR->getHcalBarrelParameters().getLayerLayout());

            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str);
            const std::string layerCodingString(m_encoder_str);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));
            const std::string staveCoding(this->GetStaveCoding(layerCodingString));

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    const edm4hep::CalorimeterHit& pCaloHit0 = pCaloHitCollection.at(i);
                    const edm4hep::CalorimeterHit* pCaloHit = &(pCaloHit0);

                    if (NULL == pCaloHit)
                        throw ("CreateHCalCaloHits Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_hCalBarrelOuterZ)
                    {
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_hCalBarrelInnerSymmetry, m_hCalBarrelInnerPhi0,
                            m_hCalBarrelInnerSymmetry - int(cellIdDecoder(pCaloHit)[ staveCoding] / 2), caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(const_cast<edm4hep::CalorimeterHit*>(pCaloHit));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout << "Error, Failed to extract hcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout<<"WARNING Failed to extract hcal calo hit: " << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract hcal calo hit collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateMuonCaloHits(const CollectionMaps& collectionMaps)
{
    for (StringVector::const_iterator iter = m_settings.m_muonCaloHitCollections.begin(), iterEnd = m_settings.m_muonCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

            std::cout<<(*iter)<<"has nElements="<<nElements<<std::endl;
            if (0 == nElements)
                continue;

            const gear::LayerLayout &endcapLayerLayout(_GEAR->getYokeEndcapParameters().getLayerLayout());
            const gear::LayerLayout &barrelLayerLayout(_GEAR->getYokeBarrelParameters().getLayerLayout()); 
            const gear::LayerLayout &plugLayerLayout(_GEAR->getYokePlugParameters().getLayerLayout());

            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str_MUON);
            const std::string layerCodingString(m_encoder_str_MUON);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));
            const std::string staveCoding(this->GetStaveCoding(layerCodingString));

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    const edm4hep::CalorimeterHit& pCaloHit0 = pCaloHitCollection.at(i);
                    const edm4hep::CalorimeterHit* pCaloHit = &(pCaloHit0);

                    if (NULL == pCaloHit)
                        throw ("Muon Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::MUON;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()] + 1;//FIXME, should use +1? it starts from 0.
                    caloHitParameters.m_isInOuterSamplingLayer = true;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    const float radius(std::sqrt(pCaloHit->getPosition()[0] * pCaloHit->getPosition()[0] +
                        pCaloHit->getPosition()[1] * pCaloHit->getPosition()[1]));

                    const bool isWithinCoil(radius < m_coilOuterR);
                    const bool isInBarrelRegion(std::fabs(pCaloHit->getPosition()[2]) < m_muonBarrelOuterZ);

                    float absorberCorrection(1.);

                    if (isInBarrelRegion && isWithinCoil)
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, plugLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else if (isInBarrelRegion)
                    {
                        this->GetBarrelCaloHitProperties(pCaloHit, barrelLayerLayout, m_muonBarrelInnerSymmetry, m_muonBarrelInnerPhi0,
                            cellIdDecoder(pCaloHit)[ staveCoding ], caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }

                    if (m_settings.m_muonDigitalHits > 0)
                    {
                        caloHitParameters.m_isDigital = true;
                        caloHitParameters.m_inputEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_hadronicEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_electromagneticEnergy = m_settings.m_muonHitEnergy;
                        caloHitParameters.m_mipEquivalentEnergy = 1.f;
                    }
                    else
                    {
                        caloHitParameters.m_isDigital = false;
                        caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_hadronicEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_electromagneticEnergy = pCaloHit->getEnergy();
                        caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_muonToMip;
                    }

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(const_cast<edm4hep::CalorimeterHit*>(pCaloHit));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout << "Failed to extract muon hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout << "Failed to extract muon hit: "  << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract muon hit collection: " << *iter  << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateLCalCaloHits(const CollectionMaps& collectionMaps)
{
    for (StringVector::const_iterator iter = m_settings.m_lCalCaloHitCollections.begin(), iterEnd = m_settings.m_lCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

            std::cout<<(*iter)<<"has nElements="<<nElements<<std::endl;
            if (0 == nElements)
                continue;

            const gear::LayerLayout &endcapLayerLayout(_GEAR->getLcalParameters().getLayerLayout()); 

            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str_LCal);
            const std::string layerCodingString(m_encoder_str_LCal);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    const edm4hep::CalorimeterHit& pCaloHit0 = pCaloHitCollection.at(i);
                    const edm4hep::CalorimeterHit* pCaloHit = &(pCaloHit0);

                    if (NULL == pCaloHit)
                        throw ("LCal Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::ECAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                        continue;

                    caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->getEnergy();
                    caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(const_cast<edm4hep::CalorimeterHit*>(pCaloHit));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout << "Failed to extract lcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout << "Failed to extract lcal calo hit: " << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract lcal calo hit collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateLHCalCaloHits(const CollectionMaps& collectionMaps)
{
    for (StringVector::const_iterator iter = m_settings.m_lHCalCaloHitCollections.begin(), iterEnd = m_settings.m_lHCalCaloHitCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

            std::cout<<(*iter)<<"has nElements="<<nElements<<std::endl;
            if (0 == nElements)
                continue;

            const gear::LayerLayout &endcapLayerLayout(_GEAR->getLHcalParameters().getLayerLayout());

            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str_LHCal);
            const std::string layerCodingString(m_encoder_str_LHCal);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));

            for (int i = 0; i < nElements; ++i)
            {
                try
                {
                    const edm4hep::CalorimeterHit& pCaloHit0 = pCaloHitCollection.at(i);
                    const edm4hep::CalorimeterHit* pCaloHit = &(pCaloHit0);

                    if (NULL == pCaloHit)
                        throw ("LHCal Collection type mismatch");

                    PandoraApi::CaloHit::Parameters caloHitParameters;
                    caloHitParameters.m_hitType = pandora::HCAL;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_layer = cellIdDecoder(pCaloHit)[layerCoding.c_str()];
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);
                    this->GetEndCapCaloHitProperties(pCaloHit, endcapLayerLayout, caloHitParameters, absorberCorrection);

                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(const_cast<edm4hep::CalorimeterHit*>(pCaloHit));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout << "Failed to extract lhcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout << "Failed to extract lhcal calo hit: " << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract lhcal calo hit collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetCommonCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
{
    const float pCaloHitPosition[3]={pCaloHit->getPosition()[0], pCaloHit->getPosition()[1], pCaloHit->getPosition()[2]};
    const pandora::CartesianVector positionVector(pCaloHitPosition[0], pCaloHitPosition[1], pCaloHitPosition[2]);

    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    caloHitParameters.m_positionVector = positionVector;
    caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
    caloHitParameters.m_pParentAddress = pCaloHit;
    caloHitParameters.m_inputEnergy = pCaloHit->getEnergy();
    caloHitParameters.m_time = pCaloHit->getTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetEndCapCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::ENDCAP;

    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() - 1));
    caloHitParameters.m_cellSize0 = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSize1 = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float radiationLength((pandora::ECAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberRadLengthECal :
        (pandora::HCAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberRadLengthHCal : m_settings.m_absorberRadLengthOther);
    const float interactionLength((pandora::ECAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberIntLengthECal :
        (pandora::HCAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberIntLengthHCal : m_settings.m_absorberIntLengthOther);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(physicalLayer));
    caloHitParameters.m_nCellRadiationLengths = radiationLength * layerAbsorberThickness;
    caloHitParameters.m_nCellInteractionLengths = interactionLength * layerAbsorberThickness;

    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        std::cout<<"WARNING CaloHitCreator::GetEndCapCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layerLayout.getNLayers(); i < iMax; ++i)
    {
        const float absorberThickness(layerLayout.getAbsorberThickness(i));

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    caloHitParameters.m_cellNormalVector = (pCaloHit->getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) :
        pandora::CartesianVector(0, 0, -1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitCreator::GetBarrelCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const gear::LayerLayout &layerLayout,
    unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
    float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::BARREL;

    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), layerLayout.getNLayers() - 1));
    caloHitParameters.m_cellSize0 = layerLayout.getCellSize0(physicalLayer);
    caloHitParameters.m_cellSize1 = layerLayout.getCellSize1(physicalLayer);
    caloHitParameters.m_cellThickness = layerLayout.getThickness(physicalLayer);

    const float radiationLength((pandora::ECAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberRadLengthECal :
        (pandora::HCAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberRadLengthHCal : m_settings.m_absorberRadLengthOther);
    const float interactionLength((pandora::ECAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberIntLengthECal :
        (pandora::HCAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberIntLengthHCal : m_settings.m_absorberIntLengthOther);

    const float layerAbsorberThickness(layerLayout.getAbsorberThickness(physicalLayer));
    caloHitParameters.m_nCellRadiationLengths = radiationLength * layerAbsorberThickness;
    caloHitParameters.m_nCellInteractionLengths = interactionLength * layerAbsorberThickness;

    //std::cout<<"m_layer="<<caloHitParameters.m_layer.Get()<<",layerLayout.getNLayers() - 1="<<layerLayout.getNLayers() - 1<<",layerAbsorberThickness="<<layerAbsorberThickness<<std::endl;
    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        std::cout<<"interactionLength="<<interactionLength<<",layerAbsorberThickness="<<layerAbsorberThickness<<",radiationLength="<<radiationLength<<",physicalLayer="<<physicalLayer<<",l1="<<caloHitParameters.m_layer.Get()<<",l2="<<layerLayout.getNLayers()-1<<std::endl;
        std::cout<<"WARNIN CaloHitCreator::GetBarrelCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layerLayout.getNLayers(); i < iMax; ++i)
    {
        const float absorberThickness(layerLayout.getAbsorberThickness(i));

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }

    if (barrelSymmetryOrder > 2)
    {
        const float phi = barrelPhi0 + (2. * M_PI * static_cast<float>(staveNumber) / static_cast<float>(barrelSymmetryOrder));
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(-std::sin(phi), std::cos(phi), 0);
    }
    else
    {
        const float pCaloHitPosition[3]={pCaloHit->getPosition()[0], pCaloHit->getPosition()[1], pCaloHit->getPosition()[2]};

        if (pCaloHitPosition[1] != 0)
        {
            const float phi = barrelPhi0 + std::atan(pCaloHitPosition[0] / pCaloHitPosition[1]);
            caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::sin(phi), std::cos(phi), 0);
        }
        else
        {
            caloHitParameters.m_cellNormalVector = (pCaloHitPosition[0] > 0) ? pandora::CartesianVector(1, 0, 0) :
                pandora::CartesianVector(-1, 0, 0);
        }
    }
}

void CaloHitCreator::GetBarrelCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer> &layerLayout,
    unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
    float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::BARREL;

    const int physicalLayer = std::max(0, int(caloHitParameters.m_layer.Get()-1) );
    caloHitParameters.m_cellSize0 = layerLayout.at(physicalLayer).cellSize0;
    caloHitParameters.m_cellSize1 = layerLayout.at(physicalLayer).cellSize1;
    caloHitParameters.m_cellThickness = layerLayout.at(physicalLayer).sensitive_thickness;

    const float radiationLength((pandora::ECAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberRadLengthECal :
        (pandora::HCAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberRadLengthHCal : m_settings.m_absorberRadLengthOther);
    const float interactionLength((pandora::ECAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberIntLengthECal :
        (pandora::HCAL == caloHitParameters.m_hitType.Get()) ? m_settings.m_absorberIntLengthHCal : m_settings.m_absorberIntLengthOther);

    const float layerAbsorberThickness(layerLayout.at(physicalLayer).absorberThickness);
    caloHitParameters.m_nCellRadiationLengths   = radiationLength   * layerAbsorberThickness;
    caloHitParameters.m_nCellInteractionLengths = interactionLength * layerAbsorberThickness;

    
    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        std::cout<<"interactionLength="<<interactionLength<<",layerAbsorberThickness="<<layerAbsorberThickness<<",radiationLength="<<radiationLength<<",physicalLayer="<<physicalLayer<<",l1="<<caloHitParameters.m_layer.Get()<<",l2="<<layerLayout.size()-1<<std::endl;
        std::cout<<"WARNIN CaloHitCreator::GetBarrelCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
    
    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layerLayout.size(); i < iMax; ++i)
    {
        const float absorberThickness(layerLayout.at(i).absorberThickness);

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

        break;
    }
    if (barrelSymmetryOrder > 2)
    {
        float phi = barrelPhi0 + (2. * M_PI * float(staveNumber) / float(barrelSymmetryOrder));
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(-std::sin(phi), std::cos(phi), 0);
        //std::cout<<"barrelPhi0="<<barrelPhi0<<",phi="<<phi<<",GetX()="<<caloHitParameters.m_cellNormalVector.Get().GetX()<<",GetY()="<<caloHitParameters.m_cellNormalVector.Get().GetY()<<",GetZ()="<<caloHitParameters.m_cellNormalVector.Get().GetZ()<<std::endl;
    }
    else
    {
        const float pCaloHitPosition[3]={pCaloHit->getPosition()[0], pCaloHit->getPosition()[1], pCaloHit->getPosition()[2]};

        if (pCaloHitPosition[1] != 0)
        {
            const float phi = barrelPhi0 + std::atan(pCaloHitPosition[0] / pCaloHitPosition[1]);
            caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::sin(phi), std::cos(phi), 0);
        }
        else
        {
            caloHitParameters.m_cellNormalVector = (pCaloHitPosition[0] > 0) ? pandora::CartesianVector(1, 0, 0) :
                pandora::CartesianVector(-1, 0, 0);
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
int CaloHitCreator::GetBarrelLayer(const edm4hep::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer> &layerLayout) const
{
    int layer = -1 ;   
    for (unsigned int i = 0, iMax = layerLayout.size(); i < iMax; ++i)
    {
        const float distance(layerLayout.at(i).distance);
        const float sensitive_thickness(layerLayout.at(i).sensitive_thickness);
        if( (distance - 0.5*sensitive_thickness) <= pCaloHit->getPosition().x && pCaloHit->getPosition().x <= (distance + 0.5*sensitive_thickness)) {layer = i+1 ; break;}
    }
    if(layer==-1)
    {
        int lmax = layerLayout.size()-1;
        std::cout<<"Error BarrelLayer, set to default 1, Hit.x="<<pCaloHit->getPosition().x<<", min_x="<<layerLayout.at(0).distance-0.5*layerLayout.at(0).sensitive_thickness<<", max_x="<<layerLayout.at(lmax).distance+0.5*layerLayout.at(lmax).sensitive_thickness<<std::endl;
        layer = 1;
    }
    return layer;
}

int CaloHitCreator::GetNLayersFromEdge(const edm4hep::CalorimeterHit *const pCaloHit) const
{
    // Calo hit coordinate calculations
    const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_hCalBarrelOuterSymmetry, m_hCalBarrelOuterPhi0));
    const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate));
    const float caloHitAbsZ(std::fabs(pCaloHit->getPosition()[2]));

    // Distance from radial outer
    float radialDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ < m_eCalBarrelOuterZ)
    {
        radialDistanceToEdge = (m_hCalBarrelOuterR - barrelMaximumRadius) / m_hCalBarrelLayerThickness;
    }
    else
    {
        radialDistanceToEdge = (m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
    }

    // Distance from rear of endcap outer
    float rearDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ >= m_eCalBarrelOuterZ)
    {
        rearDistanceToEdge = (m_hCalEndCapOuterZ - caloHitAbsZ) / m_hCalEndCapLayerThickness;
    }
    else
    {
        const float rearDistance((m_eCalBarrelOuterZ - caloHitAbsZ) / m_hCalBarrelLayerThickness);

        if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
        {
            const float overlapDistance((m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness);
            rearDistanceToEdge = std::max(rearDistance, overlapDistance);
        }
    }

    return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloHitCreator::GetMaximumRadius(const edm4hep::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
{
    
    const float pCaloHitPosition[3]={pCaloHit->getPosition()[0], pCaloHit->getPosition()[1], pCaloHit->getPosition()[2]};
    if (symmetryOrder <= 2)
        return std::sqrt((pCaloHitPosition[0] * pCaloHitPosition[0]) + (pCaloHitPosition[1] * pCaloHitPosition[1]));

    float maximumRadius(0.f);
    const float twoPi(2.f * M_PI);

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
        float radius = pCaloHitPosition[0] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

        if (radius > maximumRadius)
            maximumRadius = radius;
    }

    return maximumRadius;
}

void CaloHitCreator::GetCoding(const edm4hep::CalorimeterHit* pCaloHit, long& sys, long& x, long& y, long& z) const
{
    //sys = (pCaloHit->getCellID() << (64-8)) >> (64-8) ;
    //x   = (pCaloHit->getCellID() << (64-(8+16))) >> (64-(8+16)) ;
    //y   = (pCaloHit->getCellID() << (64-(8+16+16))) >> (64-(8+16+16)) ;
    //z   = (pCaloHit->getCellID() << (64-(8+16+16+16))) >> (64-(8+16+16+16)) ;
    unsigned long long id = pCaloHit->getCellID();
    sys =  id               &  0xFF ;
    x   = (id >> 8)         &  0xFFFF ;
    y   = (id >> (8+16))    &  0xFFFF ;
    z   = (id >> (8+16+16)) &  0xFFFF ;
}
//------------------------------------------------------------------------------------------------------------------------------------------

std::string CaloHitCreator::GetLayerCoding(const std::string &encodingString) const
{
    if (encodingString.find("layer") != std::string::npos)
        return std::string("layer");

    if (encodingString.find("K-1") != std::string::npos)
        return std::string("K-1");

    if (encodingString.find("K") != std::string::npos)
        return std::string("K");

    return std::string("unknown_layer_encoding");
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string CaloHitCreator::GetStaveCoding(const std::string &encodingString) const
{
    if (encodingString.find("stave") != std::string::npos)
        return std::string("stave");

    if (encodingString.find("S-1") != std::string::npos)
        return std::string("S-1");

    return std::string("unknown_stave_encoding") ;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitCreator::Settings::Settings() :
    m_absorberRadLengthECal(1.f),
    m_absorberIntLengthECal(1.f),
    m_absorberRadLengthHCal(1.f),
    m_absorberIntLengthHCal(1.f),
    m_absorberRadLengthOther(1.f),
    m_absorberIntLengthOther(1.f),
    m_eCalToMip(1.f),
    m_hCalToMip(1.f),
    m_muonToMip(1.f),
    m_eCalMipThreshold(0.f),
    m_hCalMipThreshold(0.f),
    m_muonMipThreshold(0.f),
    m_eCalToEMGeV(1.f),
    m_eCalToHadGeVBarrel(1.f),
    m_eCalToHadGeVEndCap(1.f),
    m_hCalToEMGeV(1.f),
    m_hCalToHadGeV(1.f),
    m_muonDigitalHits(1),
    m_muonHitEnergy(0.5f),
    m_maxHCalHitHadronicEnergy(10000.f),
    m_nOuterSamplingLayers(3),
    m_layersFromEdgeMaxRearDistance(250.f),
    m_hCalEndCapInnerSymmetryOrder(4),
    m_hCalEndCapInnerPhiCoordinate(0.f),
    m_stripSplittingOn(0),
    m_useEcalScLayers(0),
    m_eCalSiToMip(1.f),
    m_eCalScToMip(1.f),
    m_eCalSiMipThreshold(0.f),
    m_eCalScMipThreshold(0.f),
    m_eCalSiToEMGeV(1.f),
    m_eCalScToEMGeV(1.f),
    m_eCalSiToHadGeVBarrel(1.f),
    m_eCalScToHadGeVBarrel(1.f),
    m_eCalSiToHadGeVEndCap(1.f),
    m_eCalScToHadGeVEndCap(1.f)
{
}
