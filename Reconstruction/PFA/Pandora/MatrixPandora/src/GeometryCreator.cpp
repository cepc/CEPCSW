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

#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

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

//pandora::StatusCode GeometryCreator::CreateGeometry() const
//pandora::StatusCode GeometryCreator::CreateGeometry(ISvcLocator* svcloc) const
pandora::StatusCode GeometryCreator::CreateGeometry(ISvcLocator* svcloc) 
{
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
    

    //auto _gear = service<IGearSvc>("GearSvc");
    //if ( !_gear ) 
    //{
    //    throw "Failed to find GearSvc ...";
    //}
    //_GEAR = _gear->getGearMgr();
    

    try
    {
        SubDetectorTypeMap subDetectorTypeMap;
        this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

        SubDetectorNameMap subDetectorNameMap;
        this->SetAdditionalSubDetectorParameters(subDetectorNameMap);
        
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
    std::cout << "Begin SetMandatorySubDetectorParameters:" << std::endl;
    PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters, hCalBarrelParameters, hCalEndCapParameters,
        muonBarrelParameters, muonEndCapParameters;

    this->SetDefaultSubDetectorParameters(m_dd4hep, "CaloDetector", pandora::ECAL_BARREL, eCalBarrelParameters);
    //this->SetDefaultSubDetectorParameters(_GEAR->getEcalBarrelParameters(), "ECalBarrel", pandora::ECAL_BARREL, eCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(_GEAR->getEcalEndcapParameters(), "ECalEndCap", pandora::ECAL_ENDCAP, eCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(_GEAR->getHcalBarrelParameters(), "HCalBarrel", pandora::HCAL_BARREL, hCalBarrelParameters);
    this->SetDefaultSubDetectorParameters(_GEAR->getHcalEndcapParameters(), "HCalEndCap", pandora::HCAL_ENDCAP, hCalEndCapParameters);
    this->SetDefaultSubDetectorParameters(_GEAR->getYokeBarrelParameters(), "MuonBarrel", pandora::MUON_BARREL, muonBarrelParameters);
    this->SetDefaultSubDetectorParameters(_GEAR->getYokeEndcapParameters(), "MuonEndCap", pandora::MUON_ENDCAP, muonEndCapParameters);

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
    trackerParameters.m_innerRCoordinate = tpcParameters.getPadLayout().getPlaneExtent()[0];
    trackerParameters.m_innerZCoordinate = 0.f;
    trackerParameters.m_innerPhiCoordinate = 0.f;
    trackerParameters.m_innerSymmetryOrder = 0;
    trackerParameters.m_outerRCoordinate = tpcParameters.getPadLayout().getPlaneExtent()[1];
    trackerParameters.m_outerZCoordinate = tpcParameters.getMaxDriftLength();
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

void GeometryCreator::SetAdditionalSubDetectorParameters(SubDetectorNameMap &subDetectorNameMap) const
{
    std::cout << "Begin SetAdditionalSubDetectorParameters:" << std::endl;
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
         std::cout<< "warning pandora geometry creator: " << exception.what() << std::endl;
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
    std::cout << "inputParameters.getExtent()[0]="<<inputParameters.getExtent()[0]<<",inputParameters.getExtent()[1]="<<inputParameters.getExtent()[1]<<",inputParameters.getExtent()[2]="<<inputParameters.getExtent()[2]<<",inputParameters.getExtent()[3]="<<inputParameters.getExtent()[3]<< std::endl;

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
    }
}

void GeometryCreator::SetDefaultSubDetectorParameters(const dd4hep::Detector* theDetector , const std::string &subDetectorName,
    const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const
{
    //const std::vector< dd4hep::DetElement > &detElements = theDetector.detectors("EcalMatrix", true);
    const dd4hep::DetElement &detElement = theDetector->detector(subDetectorName);
    dd4hep::rec::LayeredCalorimeterData* Data = detElement.extension<dd4hep::rec::LayeredCalorimeterData>() ;
 

    parameters.m_subDetectorName    = subDetectorName;
    parameters.m_subDetectorType    = subDetectorType;
    parameters.m_innerRCoordinate   = Data->extent[0];
    parameters.m_innerZCoordinate   = Data->extent[2];
    parameters.m_innerPhiCoordinate = Data->inner_phi0;
    parameters.m_innerSymmetryOrder = Data->inner_symmetry;
    parameters.m_outerRCoordinate   = Data->extent[1];
    parameters.m_outerZCoordinate   = Data->extent[3];
    parameters.m_outerPhiCoordinate = Data->outer_phi0;
    parameters.m_outerSymmetryOrder = Data->outer_symmetry;
    parameters.m_nLayers            = Data->layers.size();
    parameters.m_isMirroredInZ      = true;
    const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>& layers= Data->layers;
    std::cout << "Extent()[0]="<<Data->extent[0]<<",Extent()[1]="<<Data->extent[1]<<",Extent()[2]="<<Data->extent[2]<<",Extent()[3]="<<Data->extent[3]<< std::endl;

    // ATTN Not always going to be correct for any optional subdetectors, but impact of this is negligible for ILD
    const float radiationLength(((pandora::ECAL_BARREL == subDetectorType) || (pandora::ECAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberRadLengthECal :
        ((pandora::HCAL_BARREL == subDetectorType) || (pandora::HCAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberRadLengthHCal : m_settings.m_absorberRadLengthOther);
    const float interactionLength(((pandora::ECAL_BARREL == subDetectorType) || (pandora::ECAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberIntLengthECal :
        ((pandora::HCAL_BARREL == subDetectorType) || (pandora::HCAL_ENDCAP == subDetectorType)) ? m_settings.m_absorberIntLengthHCal : m_settings.m_absorberIntLengthOther);
   
    //for (int i = 0; i < layerLayout.getNLayers(); ++i)
    for (auto const& layer : layers)
    {
        PandoraApi::Geometry::LayerParameters layerParameters;
        //layerParameters.m_closestDistanceToIp = layer.distance + 0.5 * layer.sensitive_thickness;//Thickness of the sensitive element (e.g. scintillator)
        layerParameters.m_closestDistanceToIp = layer.distance ;
        layerParameters.m_nRadiationLengths   = radiationLength   * layer.absorberThickness;
        layerParameters.m_nInteractionLengths = interactionLength * layer.absorberThickness;
        std::cout<<"m_closestDistanceToIp="<<layerParameters.m_closestDistanceToIp.Get()<<",m_nRadiationLengths="<<layerParameters.m_nRadiationLengths.Get()<<",m_nInteractionLengths="<<layerParameters.m_nInteractionLengths.Get()<<std::endl;
        parameters.m_layerParametersVector.push_back(layerParameters);
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
