/**
 * 
 *  @brief  Implementation of the geometry creator class.
 * 
 *  $Log: $
 */
#include "GaudiKernel/IService.h"
#include "GearSvc/IGearSvc.h"
#include "gear/BField.h"
#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/CalorimeterParameters.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/LayerLayout.h"
#include "GeometryCreator.h"
#include "Utility.h"


GeometryCreator::GeometryCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

GeometryCreator::~GeometryCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateGeometry(ISvcLocator* svcloc) 
{
    IGearSvc*  iSvc = 0;
    StatusCode sc = svcloc->service("GearSvc", iSvc, false);
    if ( !sc ) 
    {
        throw "Failed to find GearSvc ...";
    }
    _GEAR = iSvc->getGearMgr();
    

    try
    {
        SubDetectorTypeMap subDetectorTypeMap;
        this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

        SubDetectorNameMap subDetectorNameMap;
        if(m_settings.m_use_dd4hep_geo) this->SetAdditionalSubDetectorParametersDD(subDetectorNameMap);
        else                            this->SetAdditionalSubDetectorParameters  (subDetectorNameMap);
        
        if (std::string::npos != _GEAR->getDetectorName().find("ILD"))
            this->SetILDSpecificGeometry(subDetectorTypeMap, subDetectorNameMap);
        
        for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(), iterEnd = subDetectorTypeMap.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(*m_pPandora, iter->second));

        for (SubDetectorNameMap::const_iterator iter = subDetectorNameMap.begin(), iterEnd = subDetectorNameMap.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::SubDetector::Create(*m_pPandora, iter->second));
    }
    catch (gear::Exception &exception)
    {
        std::cout << "Failure in marlin pandora geometry creator, gear exception: " << exception.what() << std::endl;
        throw exception;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeometryCreator::SetMandatorySubDetectorParameters(SubDetectorTypeMap &subDetectorTypeMap) const
{
    std::cout << "GeometryCreator SetMandatorySubDetectorParameters:" << std::endl;
    PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters, hCalBarrelParameters, hCalEndCapParameters, muonBarrelParameters, muonEndCapParameters;

    if(m_settings.m_use_dd4hep_geo){
        std::cout << "Use dd4hep geo info for ECAL, HCAL and Muon" << std::endl;
        dd4hep::rec::LayeredCalorimeterData* layer_ecal_barrel = const_cast<dd4hep::rec::LayeredCalorimeterData*>(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) ) );
        if(layer_ecal_barrel) this->SetDefaultSubDetectorParameters(*layer_ecal_barrel, "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
        else{
            std::cout << "GeometryCreator:Do not find ECAL Barrel geo info from dd4hep, Try to get it from Gear" << std::endl;
            this->SetDefaultSubDetectorParameters(_GEAR->getEcalBarrelParameters(), "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
        }

        dd4hep::rec::LayeredCalorimeterData* layer_ecal_endcap = const_cast<dd4hep::rec::LayeredCalorimeterData*>(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) ) );
        if(layer_ecal_endcap) this->SetDefaultSubDetectorParameters(*layer_ecal_endcap, "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
        else{ 
            std::cout << "GeometryCreator:Do not find ECAL Endcap geo info from dd4hep, Try to get it from Gear" << std::endl;
            this->SetDefaultSubDetectorParameters(_GEAR->getEcalEndcapParameters(), "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
        }

        dd4hep::rec::LayeredCalorimeterData* layer_hcal_barrel = const_cast<dd4hep::rec::LayeredCalorimeterData*>(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) ) );
        if(layer_hcal_barrel) this->SetDefaultSubDetectorParameters(*layer_hcal_barrel, "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
        else{
            std::cout << "GeometryCreator:Do not find HCAL Barrel geo info from dd4hep, Try to get it from Gear" << std::endl;
            this->SetDefaultSubDetectorParameters(_GEAR->getHcalBarrelParameters(), "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
        }

        dd4hep::rec::LayeredCalorimeterData* layer_hcal_endcap = const_cast<dd4hep::rec::LayeredCalorimeterData*>(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) ) );
        if(layer_hcal_endcap) this->SetDefaultSubDetectorParameters(*layer_hcal_endcap, "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
        else{
            std::cout << "GeometryCreator:Do not find HCAL Endcap geo info from dd4hep, Try to get it from Gear" << std::endl;
            this->SetDefaultSubDetectorParameters(_GEAR->getHcalEndcapParameters(), "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
        }

        dd4hep::rec::LayeredCalorimeterData * layer_muon_barrel = const_cast<dd4hep::rec::LayeredCalorimeterData*>(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD ) ) );
        if(layer_muon_barrel) this->SetDefaultSubDetectorParameters(*layer_muon_barrel, "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
        else{ 
            std::cout << "GeometryCreator:Do not find Muon Barrel geo info from dd4hep, Try to get it from Gear" << std::endl;
            this->SetDefaultSubDetectorParameters(_GEAR->getYokeBarrelParameters(), "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
        }

        dd4hep::rec::LayeredCalorimeterData* layer_muon_endcap = const_cast<dd4hep::rec::LayeredCalorimeterData*>(PanUtil::getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY ) ) );
        if(layer_muon_endcap) this->SetDefaultSubDetectorParameters(*layer_muon_endcap, "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);
        else{ 
            std::cout << "GeometryCreator:Do not find Muon Endcap geo info from dd4hep, Try to get it from Gear" << std::endl;
            this->SetDefaultSubDetectorParameters(_GEAR->getYokeEndcapParameters(), "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);
        }
    }
    else{
        this->SetDefaultSubDetectorParameters(_GEAR->getEcalBarrelParameters(), "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
        this->SetDefaultSubDetectorParameters(_GEAR->getEcalEndcapParameters(), "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
        this->SetDefaultSubDetectorParameters(_GEAR->getHcalBarrelParameters(), "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
        this->SetDefaultSubDetectorParameters(_GEAR->getHcalEndcapParameters(), "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
        this->SetDefaultSubDetectorParameters(_GEAR->getYokeBarrelParameters(), "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
        this->SetDefaultSubDetectorParameters(_GEAR->getYokeEndcapParameters(), "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);
    }


    subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
    subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;
    subDetectorTypeMap[pandora::HCAL_BARREL] = hCalBarrelParameters;
    subDetectorTypeMap[pandora::HCAL_ENDCAP] = hCalEndCapParameters;
    subDetectorTypeMap[pandora::MUON_BARREL] = muonBarrelParameters;
    subDetectorTypeMap[pandora::MUON_ENDCAP] = muonEndCapParameters;

    PandoraApi::Geometry::SubDetector::Parameters trackerParameters;
    const gear::TPCParameters &tpcParameters(_GEAR->getTPCParameters());
    trackerParameters.m_subDetectorName = "Tracker";
    trackerParameters.m_subDetectorType = pandora::INNER_TRACKER;
    if(m_settings.m_use_dd4hep_geo){
        std::vector<double> tmp_extent = PanUtil::getTrackingRegionExtent();
        trackerParameters.m_innerRCoordinate = tmp_extent[0];
        trackerParameters.m_outerRCoordinate = tmp_extent[1];
        trackerParameters.m_outerZCoordinate = tmp_extent[2];
        std::cout<<"DD m_innerRCoordinate="<<trackerParameters.m_innerRCoordinate.Get()<<",m_outerRCoordinate="<<trackerParameters.m_outerRCoordinate.Get()<<",m_outerZCoordinate="<<trackerParameters.m_outerZCoordinate.Get()<<std::endl;
    }
    else{
        trackerParameters.m_innerRCoordinate = tpcParameters.getPadLayout().getPlaneExtent()[0];
        trackerParameters.m_outerRCoordinate = tpcParameters.getPadLayout().getPlaneExtent()[1];
        trackerParameters.m_outerZCoordinate = tpcParameters.getMaxDriftLength();
        std::cout<<"Gear m_innerRCoordinate="<<trackerParameters.m_innerRCoordinate.Get()<<",m_outerRCoordinate="<<trackerParameters.m_outerRCoordinate.Get()<<",m_outerZCoordinate="<<trackerParameters.m_outerZCoordinate.Get()<<std::endl;
    }
    trackerParameters.m_innerZCoordinate = 0.f;
    trackerParameters.m_innerPhiCoordinate = 0.f;
    trackerParameters.m_innerSymmetryOrder = 0;
    trackerParameters.m_outerPhiCoordinate = 0.f;
    trackerParameters.m_outerSymmetryOrder = 0;
    trackerParameters.m_isMirroredInZ = true;
    trackerParameters.m_nLayers = 0;
    subDetectorTypeMap[pandora::INNER_TRACKER] = trackerParameters;

    PandoraApi::Geometry::SubDetector::Parameters coilParameters;
    const gear::GearParameters &gearParameters(_GEAR->getGearParameters("CoilParameters"));
    coilParameters.m_subDetectorName = "Coil";
    coilParameters.m_subDetectorType = pandora::COIL;
    coilParameters.m_innerRCoordinate = gearParameters.getDoubleVal("Coil_cryostat_inner_radius");
    coilParameters.m_innerZCoordinate = 0.f;
    coilParameters.m_innerPhiCoordinate = 0.f;
    coilParameters.m_innerSymmetryOrder = 0;
    coilParameters.m_outerRCoordinate = gearParameters.getDoubleVal("Coil_cryostat_outer_radius");
    coilParameters.m_outerZCoordinate = gearParameters.getDoubleVal("Coil_cryostat_half_z");
    coilParameters.m_outerPhiCoordinate = 0.f;
    coilParameters.m_outerSymmetryOrder = 0;
    coilParameters.m_isMirroredInZ = true;
    coilParameters.m_nLayers = 0;
    subDetectorTypeMap[pandora::COIL] = coilParameters;
    std::cout << "End SetMandatorySubDetectorParameters" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
void GeometryCreator::SetAdditionalSubDetectorParametersDD(SubDetectorNameMap &subDetectorNameMap) const
{
    
  dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
  const std::vector< dd4hep::DetElement>& theECalOtherDetectors = dd4hep::DetectorSelector(mainDetector).detectors(dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::AUXILIARY ) ;
  for (std::vector< dd4hep::DetElement>::const_iterator iter = theECalOtherDetectors.begin(), iterEnd = theECalOtherDetectors.end();iter != iterEnd; ++iter){
    try
    {
        const dd4hep::DetElement& theDetector = *iter;
        const dd4hep::rec::LayeredCalorimeterData * theExtension = theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();

        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(theExtension),theDetector.name(), pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
        std::cout << "pandora geometry creator during Other ECal construction: " << exception.what() << std::endl;
    }
  }
  
  const std::vector< dd4hep::DetElement>& theHCalOtherDetectors = dd4hep::DetectorSelector(mainDetector).detectors(dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::AUXILIARY ) ;
  for (std::vector< dd4hep::DetElement>::const_iterator iter = theHCalOtherDetectors.begin(), iterEnd = theHCalOtherDetectors.end();iter != iterEnd; ++iter){
    try
    {
        const dd4hep::DetElement& theDetector = *iter;
        const dd4hep::rec::LayeredCalorimeterData * theExtension = theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();

        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(theExtension),theDetector.name(), pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
        std::cout << "pandora geometry creator during Other HCal construction: " << exception.what() << std::endl;
    }
  }
  
  const std::vector< dd4hep::DetElement>& theMuonOtherDetectors = dd4hep::DetectorSelector(mainDetector).detectors(dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON| dd4hep::DetType::AUXILIARY ) ;
  for (std::vector< dd4hep::DetElement>::const_iterator iter = theMuonOtherDetectors.begin(), iterEnd = theMuonOtherDetectors.end();iter != iterEnd; ++iter){
    try
    {
        const dd4hep::DetElement& theDetector = *iter;
        const dd4hep::rec::LayeredCalorimeterData * theExtension = theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();

        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(theExtension),theDetector.name(), pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (std::runtime_error &exception)
    {
        std::cout << "pandora geometry creator during Other Muon construction: " << exception.what() << std::endl;
    }
  }

}


void GeometryCreator::SetAdditionalSubDetectorParameters(SubDetectorNameMap &subDetectorNameMap) const
{
    std::cout << "SetAdditionalSubDetectorParameters:" << std::endl;
    try
    {
        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(_GEAR->getEcalPlugParameters(), "ECalPlug", pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (gear::Exception &exception)
    {
        std::cout << " warning pandora geometry creator: " << exception.what() << std::endl;
    }

    try
    {
        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(_GEAR->getHcalRingParameters(), "HCalRing", pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (gear::Exception &exception)
    {
         std::cout << "warning pandora geometry creator: " << exception.what() << std::endl;
    }

    try
    {
        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(_GEAR->getLcalParameters(), "LCal", pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (gear::Exception &exception)
    {
         std::cout << "warning pandora geometry creator: " << exception.what() << std::endl;
    }

    try
    {
        PandoraApi::Geometry::SubDetector::Parameters parameters;
        this->SetDefaultSubDetectorParameters(_GEAR->getLHcalParameters(), "LHCal", pandora::SUB_DETECTOR_OTHER, parameters);
        subDetectorNameMap[parameters.m_subDetectorName.Get()] = parameters;
    }
    catch (gear::Exception &exception)
    {
         std::cout << "warning pandora geometry creator: " << exception.what() << std::endl;
    }
    std::cout << "End SetAdditionalSubDetectorParameters" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
void GeometryCreator::SetDummySubDetectorParameters(const std::map<std::string, float>& para_map, const std::string &subDetectorName,
    const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
{
    parameters.m_subDetectorName = subDetectorName;
    parameters.m_subDetectorType = subDetectorType;
    parameters.m_innerRCoordinate = para_map["innerRCoordinate"];
    parameters.m_innerZCoordinate = para_map["innerZCoordinate"];
    parameters.m_innerPhiCoordinate = para_map["innerPhiCoordinate"];
    parameters.m_innerSymmetryOrder = para_map["innerSymmetryOrder"];
    parameters.m_outerRCoordinate = para_map["outerRCoordinate"];
    parameters.m_outerZCoordinate = para_map["outerZCoordinate"];
    parameters.m_outerPhiCoordinate = para_map["outerPhiCoordinate"];
    parameters.m_outerSymmetryOrder = para_map["outerSymmetryOrder"];
    parameters.m_isMirroredInZ = true;

    PandoraApi::Geometry::LayerParameters layerParameters;
    layerParameters.m_closestDistanceToIp = para_map["closestDistanceToIp"];
    layerParameters.m_nRadiationLengths = 1;
    layerParameters.m_nInteractionLengths = 1;
    parameters.m_layerParametersVector.push_back(layerParameters);
}
*/
void GeometryCreator::SetDefaultSubDetectorParameters(const gear::CalorimeterParameters &inputParameters, const std::string &subDetectorName,
    const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
{
    const gear::LayerLayout &layerLayout = inputParameters.getLayerLayout();

    parameters.m_subDetectorName = subDetectorName;
    parameters.m_subDetectorType = subDetectorType;
    parameters.m_innerRCoordinate = inputParameters.getExtent()[0];
    parameters.m_innerZCoordinate = inputParameters.getExtent()[2];
    parameters.m_innerPhiCoordinate = inputParameters.getPhi0();
    parameters.m_innerSymmetryOrder = inputParameters.getSymmetryOrder();
    parameters.m_outerRCoordinate = inputParameters.getExtent()[1];
    parameters.m_outerZCoordinate = inputParameters.getExtent()[3];
    parameters.m_outerPhiCoordinate = inputParameters.getPhi0();
    parameters.m_outerSymmetryOrder = inputParameters.getSymmetryOrder();
    parameters.m_isMirroredInZ = true;
    parameters.m_nLayers = layerLayout.getNLayers();

    std::cout <<"Gear m_subDetectorName="<<subDetectorName<<",m_subDetectorType="<<subDetectorType<<",m_innerRCoordinate="<<parameters.m_innerRCoordinate.Get()<<",m_innerZCoordinate="<<parameters.m_innerZCoordinate.Get()<<",m_innerPhiCoordinate="<<parameters.m_innerPhiCoordinate.Get()<<",m_innerSymmetryOrder="<<parameters.m_innerSymmetryOrder.Get()<<",m_outerRCoordinate="<<parameters.m_outerRCoordinate.Get()<<",m_outerZCoordinate="<<parameters.m_outerZCoordinate.Get()<<",m_outerPhiCoordinate="<<parameters.m_outerPhiCoordinate.Get()<<",m_outerSymmetryOrder="<<parameters.m_outerSymmetryOrder.Get()<<",m_nLayers="<<parameters.m_nLayers.Get()<<std::endl;
    // ATTN Not always going to be correct for any optional subdetectors, but impact of this is negligible for ILD
    const float radiationLength(((pandora::ECAL_BARREL == subDetectorType) || (pandora::ECAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberRadLengthECal :
        ((pandora::HCAL_BARREL == subDetectorType) || (pandora::HCAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberRadLengthHCal : m_settings.m_absorberRadLengthOther);
    const float interactionLength(((pandora::ECAL_BARREL == subDetectorType) || (pandora::ECAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberIntLengthECal :
        ((pandora::HCAL_BARREL == subDetectorType) || (pandora::HCAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberIntLengthHCal : m_settings.m_absorberIntLengthOther);

    for (int i = 0; i < layerLayout.getNLayers(); ++i)
    {
        PandoraApi::Geometry::LayerParameters layerParameters;
        layerParameters.m_closestDistanceToIp = layerLayout.getDistance(i) + (0.5 * (layerLayout.getThickness(i) + layerLayout.getAbsorberThickness(i)));
        layerParameters.m_nRadiationLengths = radiationLength * layerLayout.getAbsorberThickness(i);
        layerParameters.m_nInteractionLengths = interactionLength * layerLayout.getAbsorberThickness(i);
        parameters.m_layerParametersVector.push_back(layerParameters);
        std::cout <<"layer="<<i<<",m_closestDistanceToIp="<<layerParameters.m_closestDistanceToIp.Get()<<",m_nRadiationLengths="<<layerParameters.m_nRadiationLengths.Get()<<",m_nInteractionLengths="<<layerParameters.m_nInteractionLengths.Get()<<std::endl;
    }
}

void GeometryCreator::SetDefaultSubDetectorParameters(const dd4hep::rec::LayeredCalorimeterData &inputParameters, const std::string &subDetectorName,
    const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
{
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers= inputParameters.layers;
    parameters.m_subDetectorName = subDetectorName;
    parameters.m_subDetectorType = subDetectorType;
    parameters.m_innerRCoordinate = inputParameters.extent[0]/dd4hep::mm;
    parameters.m_innerZCoordinate = inputParameters.extent[2]/dd4hep::mm;
    parameters.m_innerPhiCoordinate = inputParameters.inner_phi0/dd4hep::rad;
    parameters.m_innerSymmetryOrder = inputParameters.inner_symmetry;
    parameters.m_outerRCoordinate = inputParameters.extent[1]/dd4hep::mm;
    parameters.m_outerZCoordinate = inputParameters.extent[3]/dd4hep::mm;
    parameters.m_outerPhiCoordinate = inputParameters.outer_phi0/dd4hep::rad;
    parameters.m_outerSymmetryOrder = inputParameters.outer_symmetry;
    parameters.m_isMirroredInZ = true;
    parameters.m_nLayers = layers.size();

    std::cout <<"DD m_subDetectorName="<<subDetectorName<<",m_subDetectorType="<<subDetectorType<<",m_innerRCoordinate="<<parameters.m_innerRCoordinate.Get()<<",m_innerZCoordinate="<<parameters.m_innerZCoordinate.Get()<<",m_innerPhiCoordinate="<<parameters.m_innerPhiCoordinate.Get()<<",m_innerSymmetryOrder="<<parameters.m_innerSymmetryOrder.Get()<<",m_outerRCoordinate="<<parameters.m_outerRCoordinate.Get()<<",m_outerZCoordinate="<<parameters.m_outerZCoordinate.Get()<<",m_outerPhiCoordinate="<<parameters.m_outerPhiCoordinate.Get()<<",m_outerSymmetryOrder="<<parameters.m_outerSymmetryOrder.Get()<<",m_nLayers="<<parameters.m_nLayers.Get()<<std::endl;
    for (size_t i = 0; i< layers.size(); i++)
    {
        const dd4hep::rec::LayeredCalorimeterStruct::Layer & theLayer = layers.at(i);
        
        PandoraApi::Geometry::LayerParameters layerParameters;
        
        double totalNumberOfRadLengths = theLayer.inner_nRadiationLengths;
        double totalNumberOfIntLengths = theLayer.inner_nInteractionLengths;
        
        if(i>0){
           //Add the numbers from previous layer's outer side
           totalNumberOfRadLengths += layers.at(i-1).outer_nRadiationLengths;
            totalNumberOfIntLengths += layers.at(i-1).outer_nInteractionLengths;
            
        }
        std::cout <<"layer="<<i<<",theLayer.distance="<<theLayer.distance<<std::endl; 
        layerParameters.m_closestDistanceToIp = (theLayer.distance+theLayer.inner_thickness)/dd4hep::mm; //Distance to center of sensitive element
        layerParameters.m_nRadiationLengths = totalNumberOfRadLengths;
        layerParameters.m_nInteractionLengths = totalNumberOfIntLengths;

        parameters.m_layerParametersVector.push_back(layerParameters);
        std::cout <<"layer="<<i<<",m_closestDistanceToIp="<<layerParameters.m_closestDistanceToIp.Get()<<",m_nRadiationLengths="<<layerParameters.m_nRadiationLengths.Get()<<",m_nInteractionLengths="<<layerParameters.m_nInteractionLengths.Get()<<",outer_thickness="<<theLayer.outer_thickness<<",inner_thickness="<<theLayer.inner_thickness<<",sensitive_thickness="<<theLayer.sensitive_thickness<<std::endl;
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::SetILDSpecificGeometry(SubDetectorTypeMap &subDetectorTypeMap, SubDetectorNameMap &subDetectorNameMap) const
{
    // Set positions of gaps in ILD detector and add information missing from GEAR parameters file
    try
    {
        const gear::CalorimeterParameters &hCalBarrelParameters = _GEAR->getHcalBarrelParameters();
        subDetectorTypeMap[pandora::HCAL_BARREL].m_outerPhiCoordinate = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0");
        subDetectorTypeMap[pandora::HCAL_BARREL].m_outerSymmetryOrder = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");
    }
    catch (gear::Exception &)
    {
        // aLaVideauGeometry
        return this->SetILD_SDHCALSpecificGeometry(subDetectorTypeMap);
    }

    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_eCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_eCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_eCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_eCalEndCapOuterPhiCoordinate;

    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_hCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_hCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_hCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_hCalEndCapOuterPhiCoordinate;

    subDetectorNameMap["HCalRing"].m_innerSymmetryOrder = m_settings.m_hCalRingInnerSymmetryOrder;
    subDetectorNameMap["HCalRing"].m_innerPhiCoordinate = m_settings.m_hCalRingInnerPhiCoordinate;
    subDetectorNameMap["HCalRing"].m_outerSymmetryOrder = m_settings.m_hCalRingOuterSymmetryOrder;
    subDetectorNameMap["HCalRing"].m_outerPhiCoordinate = m_settings.m_hCalRingOuterPhiCoordinate;

    // Gaps in detector active material
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelBoxGaps());
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalEndCapBoxGaps());
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateHCalBarrelConcentricGaps());

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::SetILD_SDHCALSpecificGeometry(SubDetectorTypeMap &subDetectorTypeMap) const
{
    // Non-default values (and those missing from GEAR parameters file)...
    // The following 2 parameters have no sense for Videau Geometry, set them to 0
    subDetectorTypeMap[pandora::HCAL_BARREL].m_outerPhiCoordinate = 0;
    subDetectorTypeMap[pandora::HCAL_BARREL].m_outerSymmetryOrder = 0;

    // Endcap is identical to standard ILD geometry, only HCAL barrel is different
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_eCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_eCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_eCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::ECAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_eCalEndCapOuterPhiCoordinate;

    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerSymmetryOrder = m_settings.m_hCalEndCapInnerSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_innerPhiCoordinate = m_settings.m_hCalEndCapInnerPhiCoordinate;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerSymmetryOrder = m_settings.m_hCalEndCapOuterSymmetryOrder;
    subDetectorTypeMap[pandora::HCAL_ENDCAP].m_outerPhiCoordinate = m_settings.m_hCalEndCapOuterPhiCoordinate;

    // TODO implement gaps between modules

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateHCalBarrelBoxGaps() const
{
    const std::string detectorName(_GEAR->getDetectorName());

    const gear::CalorimeterParameters &hCalBarrelParameters = _GEAR->getHcalBarrelParameters();
    const unsigned int innerSymmetryOrder(hCalBarrelParameters.getSymmetryOrder());
    const unsigned int outerSymmetryOrder(hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order"));

    if ((0 == innerSymmetryOrder) || (2 != outerSymmetryOrder / innerSymmetryOrder))
    {
        std::cout << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float innerRadius(hCalBarrelParameters.getExtent()[0]);
    const float outerRadius(hCalBarrelParameters.getExtent()[1]);
    const float outerZ(hCalBarrelParameters.getExtent()[3]);
    const float phi0(hCalBarrelParameters.getPhi0());

    const float staveGap(hCalBarrelParameters.getDoubleVal("Hcal_stave_gaps"));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, phi0, innerRadius, outerRadius,
        -outerZ, outerZ, staveGap));

    const float outerPseudoPhi0(M_PI / static_cast<float>(innerSymmetryOrder));
    const float cosOuterPseudoPhi0(std::cos(outerPseudoPhi0));

    if ((0 == outerPseudoPhi0) || (0.f == cosOuterPseudoPhi0))
    {
        std::cout << " Detector " << detectorName << " doesn't conform to expected ILD-specific geometry" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    const float middleStaveGap(hCalBarrelParameters.getDoubleVal("Hcal_middle_stave_gaps"));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(innerSymmetryOrder, outerPseudoPhi0,
        innerRadius / cosOuterPseudoPhi0, outerRadius, -outerZ, outerZ, middleStaveGap));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateHCalEndCapBoxGaps() const
{
    const gear::CalorimeterParameters &hCalEndCapParameters = _GEAR->getHcalEndcapParameters();

    const float staveGap(hCalEndCapParameters.getDoubleVal("Hcal_stave_gaps"));
    const float innerRadius(hCalEndCapParameters.getExtent()[0]);
    const float outerRadius(hCalEndCapParameters.getExtent()[1]);
    const float innerZ(hCalEndCapParameters.getExtent()[2]);
    const float outerZ(hCalEndCapParameters.getExtent()[3]);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(m_settings.m_hCalEndCapInnerSymmetryOrder,
        m_settings.m_hCalEndCapInnerPhiCoordinate, innerRadius, outerRadius, innerZ, outerZ, staveGap,
        pandora::CartesianVector(-innerRadius, 0, 0)));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateRegularBoxGaps(m_settings.m_hCalEndCapInnerSymmetryOrder,
        m_settings.m_hCalEndCapInnerPhiCoordinate, innerRadius, outerRadius, -outerZ, -innerZ, staveGap,
        pandora::CartesianVector(innerRadius, 0, 0)));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateHCalBarrelConcentricGaps() const
{
    const gear::CalorimeterParameters &hCalBarrelParameters = _GEAR->getHcalBarrelParameters();
    const float gapWidth(hCalBarrelParameters.getDoubleVal("Hcal_stave_gaps"));

    PandoraApi::Geometry::ConcentricGap::Parameters gapParameters;

    gapParameters.m_minZCoordinate = -0.5f * gapWidth;
    gapParameters.m_maxZCoordinate =  0.5f * gapWidth;
    gapParameters.m_innerRCoordinate = hCalBarrelParameters.getExtent()[0];
    gapParameters.m_innerPhiCoordinate = hCalBarrelParameters.getPhi0();
    gapParameters.m_innerSymmetryOrder = hCalBarrelParameters.getSymmetryOrder();
    gapParameters.m_outerRCoordinate = hCalBarrelParameters.getExtent()[1];
    gapParameters.m_outerPhiCoordinate = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_phi0");
    gapParameters.m_outerSymmetryOrder = hCalBarrelParameters.getIntVal("Hcal_outer_polygon_order");

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::ConcentricGap::Create(*m_pPandora, gapParameters));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreator::CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius,
    float minZ, float maxZ, float gapWidth, pandora::CartesianVector vertexOffset) const
{
    const pandora::CartesianVector basicGapVertex(pandora::CartesianVector(-0.5f * gapWidth, innerRadius, minZ) + vertexOffset);
    const pandora::CartesianVector basicSide1(gapWidth, 0, 0);
    const pandora::CartesianVector basicSide2(0, outerRadius - innerRadius, 0);
    const pandora::CartesianVector basicSide3(0, 0, maxZ - minZ);

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        const float phi = phi0 + (2. * M_PI * static_cast<float>(i) / static_cast<float>(symmetryOrder));
        const float sinPhi(std::sin(phi));
        const float cosPhi(std::cos(phi));

        PandoraApi::Geometry::BoxGap::Parameters gapParameters;

        gapParameters.m_vertex = pandora::CartesianVector(cosPhi * basicGapVertex.GetX() + sinPhi * basicGapVertex.GetY(),
            -sinPhi * basicGapVertex.GetX() + cosPhi * basicGapVertex.GetY(), basicGapVertex.GetZ());
        gapParameters.m_side1 = pandora::CartesianVector(cosPhi * basicSide1.GetX() + sinPhi * basicSide1.GetY(),
            -sinPhi * basicSide1.GetX() + cosPhi * basicSide1.GetY(), basicSide1.GetZ());
        gapParameters.m_side2 = pandora::CartesianVector(cosPhi * basicSide2.GetX() + sinPhi * basicSide2.GetY(),
            -sinPhi * basicSide2.GetX() + cosPhi * basicSide2.GetY(), basicSide2.GetZ());
        gapParameters.m_side3 = pandora::CartesianVector(cosPhi * basicSide3.GetX() + sinPhi * basicSide3.GetY(),
            -sinPhi * basicSide3.GetX() + cosPhi * basicSide3.GetY(), basicSide3.GetZ());

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::BoxGap::Create(*m_pPandora, gapParameters));
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

GeometryCreator::Settings::Settings() :
    m_absorberRadLengthECal(1.f),
    m_absorberIntLengthECal(1.f),
    m_absorberRadLengthHCal(1.f),
    m_absorberIntLengthHCal(1.f),
    m_absorberRadLengthOther(1.f),
    m_absorberIntLengthOther(1.f),
    m_eCalEndCapInnerSymmetryOrder(4),
    m_eCalEndCapInnerPhiCoordinate(0.f),
    m_eCalEndCapOuterSymmetryOrder(8),
    m_eCalEndCapOuterPhiCoordinate(0.f),
    m_hCalEndCapInnerSymmetryOrder(4),
    m_hCalEndCapInnerPhiCoordinate(0.f),
    m_hCalEndCapOuterSymmetryOrder(16),
    m_hCalEndCapOuterPhiCoordinate(0.f),
    m_hCalRingInnerSymmetryOrder(8),
    m_hCalRingInnerPhiCoordinate(0.f),
    m_hCalRingOuterSymmetryOrder(16),
    m_hCalRingOuterPhiCoordinate(0.f)
{
}
