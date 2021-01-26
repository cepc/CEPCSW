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

#include "cellIDDecoder.h"
#include "GaudiKernel/IService.h"
#include "GearSvc/IGearSvc.h"

#include "PandoraPFAlg.h"
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
    if(encoder_style==0) // LCIO style
    {
        m_encoder_str     = "M:3,S-1:3,I:9,J:9,K-1:6";
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
    if(m_settings.m_use_dd4hep_decoder){
        sc = svcloc->service("GeomSvc", m_geosvc, false);
        if (!sc) throw "Failed to find GeomSvc.";
    }
    if(m_settings.m_use_dd4hep_geo){
        const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension= PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );
        if(eCalBarrelExtension){ 
            m_eCalBarrelOuterZ        = eCalBarrelExtension->extent[3]/dd4hep::mm;
            m_eCalBarrelInnerPhi0     = eCalBarrelExtension->inner_phi0/dd4hep::rad;
            m_eCalBarrelInnerSymmetry = eCalBarrelExtension->inner_symmetry;
        }
        else{
            /*
            std::cout<<"CaloHitCreator:WARNING, Can't get ECAL geo info from dd4hep, set it to dummy value. "<<std::endl;
            m_eCalBarrelOuterZ        = 2400;
            m_eCalBarrelInnerPhi0     = 0;
            m_eCalBarrelInnerSymmetry = 8;
            */
            std::cout<<"CaloHitCreator:WARNING, Can't get ECAL geo info from dd4hep, get it from Gear. "<<std::endl;
            m_eCalBarrelOuterZ        = (_GEAR->getEcalBarrelParameters().getExtent()[3]);
            m_eCalBarrelInnerPhi0     = (_GEAR->getEcalBarrelParameters().getPhi0());
            m_eCalBarrelInnerSymmetry = (_GEAR->getEcalBarrelParameters().getSymmetryOrder());
        }
        //Get HCal Barrel extension by type, ignore plugs and rings 
        const dd4hep::rec::LayeredCalorimeterData * hCalBarrelExtension= PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) );
        //Get HCal Endcap extension by type, ignore plugs and rings 
        const dd4hep::rec::LayeredCalorimeterData * hCalEndcapExtension= PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) );
        if(hCalBarrelExtension){
            m_hCalBarrelOuterZ             =   hCalBarrelExtension->extent[3]/dd4hep::mm;
            m_hCalBarrelInnerPhi0          =   hCalBarrelExtension->inner_phi0/dd4hep::rad;
            m_hCalBarrelInnerSymmetry      =   hCalBarrelExtension->inner_symmetry;
            m_hCalBarrelOuterR             =   hCalBarrelExtension->extent[1]/dd4hep::mm;
            m_hCalBarrelOuterPhi0          =   hCalBarrelExtension->outer_phi0/dd4hep::rad;
            m_hCalBarrelOuterSymmetry      =   hCalBarrelExtension->outer_symmetry;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers= hCalBarrelExtension->layers;
            //Take thicknesses from last layer (was like that before with gear)
            m_hCalBarrelLayerThickness =(barrelLayers.back().inner_thickness+barrelLayers.back().outer_thickness)/dd4hep::mm;
        }
        else{
            /*
            std::cout<<"CaloHitCreator:WARNING, Can't get HCAL Barrel geo info from dd4hep, set it to dummy value. "<<std::endl;
            m_hCalBarrelOuterZ             =   0;
            m_hCalBarrelInnerPhi0          =   0;
            m_hCalBarrelInnerSymmetry      =   0;
            m_hCalBarrelOuterR             =   0;
            m_hCalBarrelOuterPhi0          =   0;
            m_hCalBarrelOuterSymmetry      =   0;
            m_hCalBarrelLayerThickness     =   1;
            */
            std::cout<<"CaloHitCreator:WARNING, Can't get HCAL Barrel geo info from dd4hep, get it from Gear"<<std::endl;
            m_hCalBarrelOuterZ        = (_GEAR->getHcalBarrelParameters().getExtent()[3]);
            m_hCalBarrelInnerPhi0     = (_GEAR->getHcalBarrelParameters().getPhi0());
            m_hCalBarrelInnerSymmetry = (_GEAR->getHcalBarrelParameters().getSymmetryOrder());
            m_hCalBarrelOuterR        = (_GEAR->getHcalBarrelParameters().getExtent()[1]);
            m_hCalBarrelOuterPhi0     = ((std::find(_GEAR->getHcalBarrelParameters().getIntKeys().begin(),
                _GEAR->getHcalBarrelParameters().getIntKeys().end(),
                "Hcal_outer_polygon_phi0") != _GEAR->getHcalBarrelParameters().getIntKeys().end() ?
                _GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0")
                : 0));
            m_hCalBarrelOuterSymmetry = ((std::find(_GEAR->getHcalBarrelParameters().getIntKeys().begin(),
                _GEAR->getHcalBarrelParameters().getIntKeys().end(),
                "Hcal_outer_polygon_order") != _GEAR->getHcalBarrelParameters().getIntKeys().end() ?
                _GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order")
                : 0));
            const gear::LayerLayout &hCalBarrelLayerLayout(_GEAR->getHcalBarrelParameters().getLayerLayout()); 
            m_hCalBarrelLayerThickness = hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1);
        }
        if(hCalEndcapExtension){
            m_hCalEndCapOuterR             =   hCalEndcapExtension->extent[1]/dd4hep::mm;
            m_hCalEndCapOuterZ             =   hCalEndcapExtension->extent[3]/dd4hep::mm;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers= hCalEndcapExtension->layers;
            //Take thicknesses from last layer (was like that before with gear)
            m_hCalEndCapLayerThickness =(endcapLayers.back().inner_thickness+endcapLayers.back().outer_thickness)/dd4hep::mm;
        }
        else{
            /*
            std::cout<<"CaloHitCreator:WARNING, Can't get HCAL Endcap geo info from dd4hep, set it to dummy value. "<<std::endl;
            m_hCalEndCapOuterR             =   0;
            m_hCalEndCapOuterZ             =   0;
            m_hCalEndCapLayerThickness     =   1;
            */
            std::cout<<"CaloHitCreator:WARNING, Can't get HCAL Endcap geo info from dd4hep, get it from Gear. "<<std::endl;
            m_hCalEndCapOuterR        = (_GEAR->getHcalEndcapParameters().getExtent()[1]);
            m_hCalEndCapOuterZ        = (_GEAR->getHcalEndcapParameters().getExtent()[3]);
            const gear::LayerLayout &hCalEndCapLayerLayout(_GEAR->getHcalEndcapParameters().getLayerLayout());
            m_hCalEndCapLayerThickness = hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1);

        }
        //Get Muon Barrel extension by type, ignore plugs and rings 
        const dd4hep::rec::LayeredCalorimeterData * muonBarrelExtension= PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) );
        //fg: muon endcap is not used :
        //Get Muon Endcap extension by type, ignore plugs and rings 
        // const dd4hep::rec::LayeredCalorimeterData * muonEndcapExtension= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY ) );
        if(muonBarrelExtension){
            m_muonBarrelOuterZ             =   muonBarrelExtension->extent[3]/dd4hep::mm;
            m_muonBarrelInnerPhi0          =   muonBarrelExtension->inner_phi0/dd4hep::rad;
            m_muonBarrelInnerSymmetry      =   muonBarrelExtension->inner_symmetry;
        }
        else{
            /*
            std::cout<<"CaloHitCreator:WARNING, Can't get Muon Barrel geo info from dd4hep, set it to dummy value. "<<std::endl;
            m_muonBarrelOuterZ             =   0;
            m_muonBarrelInnerPhi0          =   0;
            m_muonBarrelInnerSymmetry      =   0;
            */
            std::cout<<"CaloHitCreator:WARNING, Can't get Muon Barrel geo info from dd4hep, get it from Gear. "<<std::endl;
            m_muonBarrelInnerPhi0     = (_GEAR->getYokeBarrelParameters().getPhi0());
            m_muonBarrelOuterZ        = (_GEAR->getYokeBarrelParameters().getExtent()[3]);
            m_muonBarrelInnerSymmetry = (_GEAR->getYokeBarrelParameters().getSymmetryOrder());
        }
        //Get COIL extension
        const dd4hep::rec::LayeredCalorimeterData * coilExtension= PanUtil::getExtension( ( dd4hep::DetType::COIL ) );
        if(coilExtension){
            m_coilOuterR                   =   coilExtension->extent[1]/dd4hep::mm;
        }
        else{
            /*
            std::cout<<"CaloHitCreator:WARNING, Can't get Coil geo info from dd4hep, set it to dummy value. "<<std::endl;
            m_coilOuterR = 0;
            */
            std::cout<<"CaloHitCreator:WARNING, Can't get Coil geo info from dd4hep, get it from Gear. "<<std::endl;
            m_coilOuterR              = (_GEAR->getGearParameters("CoilParameters").getDoubleVal("Coil_cryostat_outer_radius"));
        }
    }
    else{
        m_eCalBarrelOuterZ        = (_GEAR->getEcalBarrelParameters().getExtent()[3]);
        m_eCalBarrelInnerPhi0     = (_GEAR->getEcalBarrelParameters().getPhi0());
        m_eCalBarrelInnerSymmetry = (_GEAR->getEcalBarrelParameters().getSymmetryOrder());
        m_hCalBarrelOuterZ        = (_GEAR->getHcalBarrelParameters().getExtent()[3]);
        m_hCalBarrelInnerPhi0     = (_GEAR->getHcalBarrelParameters().getPhi0());
        m_hCalBarrelInnerSymmetry = (_GEAR->getHcalBarrelParameters().getSymmetryOrder());
        m_hCalBarrelOuterR        = (_GEAR->getHcalBarrelParameters().getExtent()[1]);
        m_hCalBarrelOuterPhi0     = ((std::find(_GEAR->getHcalBarrelParameters().getIntKeys().begin(),
            _GEAR->getHcalBarrelParameters().getIntKeys().end(),
            "Hcal_outer_polygon_phi0") != _GEAR->getHcalBarrelParameters().getIntKeys().end() ?
            _GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0")
            : 0));
        m_hCalBarrelOuterSymmetry = ((std::find(_GEAR->getHcalBarrelParameters().getIntKeys().begin(),
            _GEAR->getHcalBarrelParameters().getIntKeys().end(),
            "Hcal_outer_polygon_order") != _GEAR->getHcalBarrelParameters().getIntKeys().end() ?
            _GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order")
            : 0));
        const gear::LayerLayout &hCalBarrelLayerLayout(_GEAR->getHcalBarrelParameters().getLayerLayout()); 
        m_hCalBarrelLayerThickness = hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1);
        m_hCalEndCapOuterR        = (_GEAR->getHcalEndcapParameters().getExtent()[1]);
        m_hCalEndCapOuterZ        = (_GEAR->getHcalEndcapParameters().getExtent()[3]);
        const gear::LayerLayout &hCalEndCapLayerLayout(_GEAR->getHcalEndcapParameters().getLayerLayout());
        m_hCalEndCapLayerThickness = hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1);
        m_coilOuterR              = (_GEAR->getGearParameters("CoilParameters").getDoubleVal("Coil_cryostat_outer_radius"));
        m_muonBarrelInnerPhi0     = (_GEAR->getYokeBarrelParameters().getPhi0());
        m_muonBarrelOuterZ        = (_GEAR->getYokeBarrelParameters().getExtent()[3]);
        m_muonBarrelInnerSymmetry = (_GEAR->getYokeBarrelParameters().getSymmetryOrder());
    }

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
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalCaloHits (collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMuonCaloHits (collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLCalCaloHits (collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateLHCalCaloHits(collectionMaps));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateECalCaloHits(const CollectionMaps& collectionMaps)
{
    for (unsigned int k = 0; k < m_settings.m_eCalCaloHitCollections.size(); k++)
    {
        std::string tmp_col_name = m_settings.m_eCalCaloHitCollections.at(k);
        if(collectionMaps.collectionMap_CaloHit.find(tmp_col_name) == collectionMaps.collectionMap_CaloHit.end()) { if(m_settings.m_debug) std::cout<<"not find "<<tmp_col_name<<std::endl; continue;}
        try
        {
            if(m_settings.m_debug) std::cout<<"CaloHitCreator for "<<tmp_col_name<<std::endl;
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(tmp_col_name))->second;
            const int nElements(pCaloHitCollection.size());

            if (0 == nElements)
                continue;


            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str);
            const std::string layerCodingString(m_encoder_str);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));
            const std::string staveCoding(this->GetStaveCoding(layerCodingString));
            // get the DD4hep readout
            const std::string name_readout = m_settings.m_eCalCaloReadOuts.at(k);
            if(m_settings.m_debug) std::cout<<"readout= "<<name_readout<<std::endl;
            if( m_settings.m_use_dd4hep_decoder ){
                m_decoder = m_geosvc->getDecoder(name_readout);
                if (!m_decoder) throw "Failed to get the decoder. ";
            }
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* barrelLayers=nullptr;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* endcapLayers=nullptr;
            const gear::LayerLayout* barrelLayerLayout=nullptr;
            const gear::LayerLayout* endcapLayerLayout=nullptr;
            if(m_settings.m_use_dd4hep_geo){
                barrelLayers= &(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )->layers);
                endcapLayers= &(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) )->layers);
            }
            else{
                barrelLayerLayout = &(_GEAR->getEcalBarrelParameters().getLayerLayout()); 
                endcapLayerLayout = &(_GEAR->getEcalEndcapParameters().getLayerLayout());
            }

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
                        std::string collectionName(tmp_col_name);
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
                    caloHitParameters.m_layer = m_settings.m_use_dd4hep_decoder == false ? cellIdDecoder(pCaloHit)[layerCoding.c_str()] + 1 : m_decoder->get(pCaloHit->getCellID(), "layer");// from 0 to 29, 0 is preshower layer
                    int Stave = 0 ; 
                    if (m_settings.m_use_dd4hep_decoder == false){
                        Stave = cellIdDecoder(pCaloHit)[ staveCoding];
                    }
                    else{
                        Stave = m_decoder->get(pCaloHit->getCellID(), "stave");
                        Stave = Stave <=2 ? Stave+5 : Stave-3 ;//change to correct style
                    }
                    //std::cout<<"0Stave="<<Stave<<",0layer="<<caloHitParameters.m_layer.Get()<<std::endl;
                    if (Stave<0) throw "wrong Stave";
                    if (m_settings.m_use_preshower==false && caloHitParameters.m_layer.Get()<1) continue;//don't use preshower layer 
                    //std::cout<<"Stave="<<Stave<<",layer="<<caloHitParameters.m_layer.Get()<<std::endl;
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_eCalBarrelOuterZ)
                    {
                        if(m_settings.m_use_dd4hep_geo) this->GetBarrelCaloHitProperties(pCaloHit, *barrelLayers, m_eCalBarrelInnerSymmetry, m_eCalBarrelInnerPhi0, Stave, caloHitParameters, absorberCorrection);
                        else                            this->GetBarrelCaloHitProperties(pCaloHit, *barrelLayerLayout, m_eCalBarrelInnerSymmetry, m_eCalBarrelInnerPhi0, Stave, caloHitParameters, absorberCorrection);
                        caloHitParameters.m_hadronicEnergy = eCalToHadGeVBarrel * pCaloHit->getEnergy();
                    }
                    else
                    {
                        if(m_settings.m_use_dd4hep_geo) this->GetEndCapCaloHitProperties(pCaloHit, *endcapLayers, caloHitParameters, absorberCorrection);
                        
                        else                            this->GetEndCapCaloHitProperties(pCaloHit, *endcapLayerLayout, caloHitParameters, absorberCorrection);
                        caloHitParameters.m_hadronicEnergy = eCalToHadGeVEndCap * pCaloHit->getEnergy();
                    }

                    //caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * eCalToMip * absorberCorrection;
                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * eCalToMip;//is absorberCorrection it needed for digi input, also the m_mipEquivalentEnergy seems is not used in alg

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
                    std::cout<<"Failed to extract ecal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout<<"Failed to extract ecal calo hit" <<  std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout<< "Failed to extract ecal calo hit collection: " << tmp_col_name << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateHCalCaloHits(const CollectionMaps& collectionMaps)
{
    for (unsigned int k = 0; k < m_settings.m_hCalCaloHitCollections.size(); k++)
    {
        std::string tmp_col_name = m_settings.m_hCalCaloHitCollections.at(k);
        if(collectionMaps.collectionMap_CaloHit.find(tmp_col_name) == collectionMaps.collectionMap_CaloHit.end()) { if(m_settings.m_debug) std::cout<<"not find "<<tmp_col_name<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(tmp_col_name))->second;
            const int nElements(pCaloHitCollection.size());

            if (0 == nElements)
                continue;


            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str);
            const std::string layerCodingString(m_encoder_str);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));
            const std::string staveCoding(this->GetStaveCoding(layerCodingString));
            // get the DD4hep readout
            const std::string name_readout = m_settings.m_hCalCaloReadOuts.at(k);
            if(m_settings.m_debug) std::cout<<"readout= "<<name_readout<<std::endl;
            if( m_settings.m_use_dd4hep_decoder ){
                m_decoder = m_geosvc->getDecoder(name_readout);
                if (!m_decoder) throw "Failed to get the decoder. ";
            }
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* barrelLayers=nullptr;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* endcapLayers=nullptr;
            const gear::LayerLayout* barrelLayerLayout=nullptr;
            const gear::LayerLayout* endcapLayerLayout=nullptr;
            if(m_settings.m_use_dd4hep_geo){
                barrelLayers= &(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) )->layers);
                endcapLayers= &(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) )->layers);
            }
            else{
                barrelLayerLayout = &(_GEAR->getHcalBarrelParameters().getLayerLayout());
                endcapLayerLayout = &(_GEAR->getHcalEndcapParameters().getLayerLayout());
            }

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
                    caloHitParameters.m_isDigital = false;// if it is DHCAL or AHCAL
                    caloHitParameters.m_layer = m_settings.m_use_dd4hep_decoder == false ? cellIdDecoder(pCaloHit)[layerCoding.c_str()] : m_decoder->get(pCaloHit->getCellID(), "layer");
                    caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);
                    int Stave = 0 ; 
                    if (m_settings.m_use_dd4hep_decoder == false){
                        Stave = cellIdDecoder(pCaloHit)[ staveCoding];
                    }
                    else{
                        Stave = m_decoder->get(pCaloHit->getCellID(), "stave");
                        Stave = Stave <=2 ? Stave+5 : Stave-3 ;//FIXME , need check!!
                    }

                    float absorberCorrection(1.);

                    if (std::fabs(pCaloHit->getPosition()[2]) < m_hCalBarrelOuterZ)
                    {
                        if(m_settings.m_use_dd4hep_geo) this->GetBarrelCaloHitProperties(pCaloHit, *barrelLayers, m_hCalBarrelInnerSymmetry, m_hCalBarrelInnerPhi0, m_hCalBarrelInnerSymmetry - int(Stave / 2), caloHitParameters, absorberCorrection);
                        else                            this->GetBarrelCaloHitProperties(pCaloHit, *barrelLayerLayout, m_hCalBarrelInnerSymmetry, m_hCalBarrelInnerPhi0, m_hCalBarrelInnerSymmetry - int(Stave / 2), caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        if(m_settings.m_use_dd4hep_geo)   this->GetEndCapCaloHitProperties(pCaloHit, *endcapLayers, caloHitParameters, absorberCorrection);
                        else                              this->GetEndCapCaloHitProperties(pCaloHit, *endcapLayerLayout, caloHitParameters, absorberCorrection);
                    }

                    //caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;
                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip;

                    if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
                        continue;

                    caloHitParameters.m_hadronicEnergy = std::min(m_settings.m_hCalToHadGeV * pCaloHit->getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
                    caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * pCaloHit->getEnergy();

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
                    m_calorimeterHitVector.push_back(const_cast<edm4hep::CalorimeterHit*>(pCaloHit));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout << "Failed to extract hcal calo hit: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout<<"Failed to extract hcal calo hit" << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract hcal calo hit collection: " << tmp_col_name << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CaloHitCreator::CreateMuonCaloHits(const CollectionMaps& collectionMaps)
{
    for (unsigned int k = 0; k < m_settings.m_muonCaloHitCollections.size(); k++)
    {
        std::string tmp_col_name = m_settings.m_muonCaloHitCollections.at(k);
        if(collectionMaps.collectionMap_CaloHit.find(tmp_col_name) == collectionMaps.collectionMap_CaloHit.end()) {if(m_settings.m_debug) std::cout<<"not find "<<tmp_col_name<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(tmp_col_name))->second;
            const int nElements(pCaloHitCollection.size());

            if (0 == nElements)
                continue;


            ID_UTIL::CellIDDecoder<const edm4hep::CalorimeterHit> cellIdDecoder(m_encoder_str_MUON);
            const std::string layerCodingString(m_encoder_str_MUON);
            const std::string layerCoding(this->GetLayerCoding(layerCodingString));
            const std::string staveCoding(this->GetStaveCoding(layerCodingString));
            // get the DD4hep readout
            const std::string name_readout = m_settings.m_muonCalCaloReadOuts.at(k);
            if(m_settings.m_debug) std::cout<<"readout= "<<name_readout<<std::endl;
            if( m_settings.m_use_dd4hep_decoder ){
                m_decoder = m_geosvc->getDecoder(name_readout);
                if (!m_decoder) throw "Failed to get the decoder. ";
            }
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* barrelLayers=nullptr;
            const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* endcapLayers=nullptr;
            const gear::LayerLayout* barrelLayerLayout=nullptr;
            const gear::LayerLayout* plugLayerLayout=nullptr;
            const gear::LayerLayout* endcapLayerLayout=nullptr;
            if(m_settings.m_use_dd4hep_geo){
                barrelLayers= &(PanUtil::getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON| dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD  ))->layers);
                endcapLayers = &(PanUtil::getExtension(( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON| dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD  ))->layers);
            }
            else{
                plugLayerLayout= &(_GEAR->getYokePlugParameters().getLayerLayout());
                barrelLayerLayout = &(_GEAR->getYokeBarrelParameters().getLayerLayout()); 
                endcapLayerLayout = &(_GEAR->getYokeEndcapParameters().getLayerLayout());
            }

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
                    caloHitParameters.m_layer = m_settings.m_use_dd4hep_decoder == false ? cellIdDecoder(pCaloHit)[layerCoding.c_str()] + 1 : m_decoder->get(pCaloHit->getCellID(), "layer");
                    caloHitParameters.m_isInOuterSamplingLayer = true;
                    this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);
                    int Stave = 0 ; 
                    if (m_settings.m_use_dd4hep_decoder == false){
                        Stave = cellIdDecoder(pCaloHit)[ staveCoding];
                    }
                    else{
                        Stave = m_decoder->get(pCaloHit->getCellID(), "stave");
                        Stave = Stave <=2 ? Stave+5 : Stave-3 ;//FIXME , need check!!
                    }

                    const float radius(std::sqrt(pCaloHit->getPosition()[0] * pCaloHit->getPosition()[0] +
                        pCaloHit->getPosition()[1] * pCaloHit->getPosition()[1]));

                    const bool isWithinCoil(radius < m_coilOuterR);
                    const bool isInBarrelRegion(std::fabs(pCaloHit->getPosition()[2]) < m_muonBarrelOuterZ);

                    float absorberCorrection(1.);

                    if (isInBarrelRegion && isWithinCoil)
                    {
                        if(m_settings.m_use_dd4hep_geo) std::cout<<"BIG WARNING: CANNOT HANDLE PLUG HITS (no plug), DO NOTHING!"<<std::endl;
                        else                            this->GetEndCapCaloHitProperties(pCaloHit, *plugLayerLayout, caloHitParameters, absorberCorrection);
                    }
                    else if (isInBarrelRegion)
                    {   if(m_settings.m_use_dd4hep_geo) this->GetBarrelCaloHitProperties(pCaloHit, *barrelLayers     , m_muonBarrelInnerSymmetry, m_muonBarrelInnerPhi0, Stave, caloHitParameters, absorberCorrection);
                        else                            this->GetBarrelCaloHitProperties(pCaloHit, *barrelLayerLayout, m_muonBarrelInnerSymmetry, m_muonBarrelInnerPhi0, Stave, caloHitParameters, absorberCorrection);
                    }
                    else
                    {
                        if(m_settings.m_use_dd4hep_geo) this->GetEndCapCaloHitProperties(pCaloHit, *endcapLayers     , caloHitParameters, absorberCorrection);
                        else                            this->GetEndCapCaloHitProperties(pCaloHit, *endcapLayerLayout, caloHitParameters, absorberCorrection);
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
                    std::cout << "Failed to extract muon hit"  << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract muon hit collection: " << tmp_col_name  << std::endl;
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
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) {if(m_settings.m_debug)std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

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

                    //caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip * absorberCorrection;
                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_eCalToMip;

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
                    std::cout << "Failed to extract lcal calo hit" << std::endl;
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
        if(collectionMaps.collectionMap_CaloHit.find(*iter) == collectionMaps.collectionMap_CaloHit.end()) {if(m_settings.m_debug) std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            const std::vector<edm4hep::CalorimeterHit>& pCaloHitCollection = (collectionMaps.collectionMap_CaloHit.find(*iter))->second;
            const int nElements(pCaloHitCollection.size());

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

                    //caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip * absorberCorrection;
                    caloHitParameters.m_mipEquivalentEnergy = pCaloHit->getEnergy() * m_settings.m_hCalToMip;

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
                    std::cout << "Failed to extract lhcal calo hit" << std::endl;
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
void CaloHitCreator::GetEndCapCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
    PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::ENDCAP;

    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
    caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0/dd4hep::mm;
    caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1/dd4hep::mm;
    //std::cout<<"Endcap m_cellSize0="<<caloHitParameters.m_cellSize0.Get()<<",m_cellSize1="<<caloHitParameters.m_cellSize1.Get()<<std::endl;
    double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
    double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
    double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;

    if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;

    }
    caloHitParameters.m_cellThickness = thickness;
    caloHitParameters.m_nCellRadiationLengths = nRadLengths;
    caloHitParameters.m_nCellInteractionLengths = nIntLengths;
    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        std::cout<<"WARNING CaloHitCreator::GetEndCapCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    absorberCorrection = 1.;
    for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
    {
        float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 )/dd4hep::mm);
        
        if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0)/dd4hep::mm;

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

    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
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

void CaloHitCreator::GetBarrelCaloHitProperties(const edm4hep::CalorimeterHit *const pCaloHit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &layers,
    unsigned int barrelSymmetryOrder, float barrelPhi0, unsigned int staveNumber, PandoraApi::CaloHit::Parameters &caloHitParameters,
    float &absorberCorrection) const
{
    caloHitParameters.m_hitRegion = pandora::BARREL;
    const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
    caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0/dd4hep::mm;
    caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1/dd4hep::mm;
    if(m_settings.m_debug) std::cout<<"DD m_cellSize0="<<caloHitParameters.m_cellSize0.Get()<<",m_cellSize1="<<caloHitParameters.m_cellSize1.Get()<<std::endl;
    double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
    double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

    double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0)/dd4hep::mm;
    }
    
    caloHitParameters.m_cellThickness = thickness;
    caloHitParameters.m_nCellRadiationLengths = nRadLengths;
    caloHitParameters.m_nCellInteractionLengths = nIntLengths;

    if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
    {
        std::cout<<"WARNIN CaloHitCreator::GetBarrelCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
    absorberCorrection = 1.;
    float absorberThickness_0 = 0; 
    for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
    {
        float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 )/dd4hep::mm);
        
        if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0)/dd4hep::mm;

        if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

        if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;
        absorberThickness_0 = absorberThickness;
        break;
    }


    //std::cout<<"DD m_cellSize0="<<caloHitParameters.m_cellSize0.Get()<<",m_cellSize1="<<caloHitParameters.m_cellSize1.Get()<<",physicalLayer="<<physicalLayer<<",layerAbsorberThickness="<<layerAbsorberThickness<<",absorberThickness_0="<<absorberThickness_0<<",m_cellThickness="<<caloHitParameters.m_cellThickness.Get()<<",m_nCellRadiationLengths="<<caloHitParameters.m_nCellRadiationLengths.Get()<<",m_nCellInteractionLengths="<<caloHitParameters.m_nCellInteractionLengths.Get()<<",absorberCorrection="<<absorberCorrection<<std::endl;


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

//------------------------------------------------------------------------------------------------------------------------------------------

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
