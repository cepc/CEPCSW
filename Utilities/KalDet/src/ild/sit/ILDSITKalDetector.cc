
#include "ILDSITKalDetector.h"

#include "MaterialDataBase.h"

#include "ILDParallelPlanarStripMeasLayer.h"


#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"

#include <gear/GEAR.h>
#include "gear/BField.h"
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>
#include "gearimpl/Util.h"

#include "TMath.h"

#include "math.h"
#include <sstream>

// #include "streamlog/streamlog.h"

ILDSITKalDetector::ILDSITKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc )
: TVKalDetector(300) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  // std::cout << "ILDSITKalDetector building SIT detector using GEAR " << std::endl ;
  
  MaterialDataBase::Instance().registerForService(gearMgr, geoSvc);
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon");

  if(geoSvc){
    this->setupGearGeom(geoSvc);
  }
  else{
    this->setupGearGeom(gearMgr) ;
  }
  
  if (_isStripDetector) {
    // streamlog_out(DEBUG4) << "\t\t building SIT detector as STRIP Detector." << std::endl ;    
  } else {
    // streamlog_out(DEBUG4) << "\t\t building SIT detector as PIXEL Detector." << std::endl ;    
  }
  
  //--The Ladder structure (realistic ladder)--
  int nLadders;
  
  //  Bool_t active = true;
  Bool_t dummy  = false;
  
  static const double eps_layer  = 1e-6; 
  static const double eps_sensor = 1e-8; 
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for (int layer=0; layer<_nLayers; ++layer) {
    
    double offset = _SITgeo[layer].offset ;
    
    if( offset!=0.0 ) {
      std::cerr << "SIT can not have offsets using the current implementation: exit(1)" << std::endl ; 
      exit(1) ;
    }
    
    nLadders = _SITgeo[layer].nLadders ;
    
    const double phi0 = _SITgeo[layer].phi0 ;
    
    const double ladder_distance = _SITgeo[layer].supRMin ;
    const double ladder_thickness = _SITgeo[layer].supThickness ;
    
    const double sensitive_distance = _SITgeo[layer].senRMin ;
    const double sensitive_thickness = _SITgeo[layer].senThickness ;
    
    const double width = _SITgeo[layer].width ;
    const double length = _SITgeo[layer].length;
    
    double currPhi;
    const double dphi = _SITgeo[layer].dphi ;
    
    const double stripAngle = pow(-1,layer) *_SITgeo[layer].stripAngle;
    
    const int nsensors = _SITgeo[layer].nSensorsPerLadder;
    
    const double sensor_length = _SITgeo[layer].sensorLength;
    
    
    for (int ladder=0; ladder<nLadders; ++ladder) {
      
      
      currPhi = phi0 + (dphi * ladder);
      
      encoder.reset() ;  // reset to 0
      
      encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::SIT ;
      encoder[lcio::ILDCellID0::side] = 0 ;
      encoder[lcio::ILDCellID0::layer]  = layer ;
      encoder[lcio::ILDCellID0::module] = ladder ;
      
      double z_centre_support = 0.0;
      
      // check if the sensitive is inside or outside for the support 
      if( sensitive_distance < ladder_distance  ) {
        
        double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+0) * eps_layer ;
        double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+2) * eps_layer ;
        double sup_back_sorting_policy          = ladder_distance     + (4 * ladder+3) * eps_layer ;
        
        // air - sensitive boundary
        Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SITSenFront")) ;
        
        for (int isensor=0; isensor<nsensors; ++isensor) {

          encoder[lcio::ILDCellID0::sensor] = isensor ;          
          int CellID = encoder.lowWord() ;
          
          double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+1) * eps_layer + eps_sensor * isensor ;
          
          double z_centre_sensor = -0.5*length + (0.5*sensor_length) + (isensor*sensor_length) ;

          
          if (_isStripDetector) {
          
            // measurement plane defined as the middle of the sensitive volume 
            Add(new ILDParallelPlanarStripMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, stripAngle, CellID, "SITMeaslayer" )) ;
            
          } else {
            // measurement plane defined as the middle of the sensitive volume
            Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, true, CellID, "SITMeaslayer" )) ;
          }
          
          
	  //std::cout << "ILDSITKalDetector add surface with CellID = " << CellID << std::endl ;
          
        }
        
       
        // sensitive - support boundary 
        Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SITSenSupportIntf" )) ; 
        
        // support - air boundary
        Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, sup_back_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SITSupRear" )) ; 
      }
      else {
        
        double sup_front_sorting_policy         = ladder_distance     + (4 * ladder+0) * eps_layer ;
        double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+1) * eps_layer ;
        double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+3) * eps_layer ;
        
        // air - support boundary
        Add(new ILDParallelPlanarMeasLayer(air, carbon, ladder_distance, currPhi, _bZ, sup_front_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SITSupFront")) ;
        
        // support boundary - sensitive
        Add(new ILDParallelPlanarMeasLayer(carbon, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SITSenSupportIntf" )) ; 
        
        
        for (int isensor=0; isensor<nsensors; ++isensor) {

          encoder[lcio::ILDCellID0::sensor] = isensor ;          
          int CellID = encoder.lowWord() ;
          
          double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+2) * eps_layer + eps_sensor * isensor ;

          double z_centre_sensor = -0.5*length + (0.5*sensor_length) + (isensor*sensor_length) ;

          
          if (_isStripDetector) {
            
          // measurement plane defined as the middle of the sensitive volume 
          Add(new ILDParallelPlanarStripMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, stripAngle, CellID, "SITMeaslayer" )) ;
            
          } else {
            // measurement plane defined as the middle of the sensitive volume
            Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, true, CellID, "SITMeaslayer" )) ;
          }

	  
	  //std::cout << "ILDSITKalDetector add surface with CellID = " << CellID << std::endl ;

          
        }

                
        // support - air boundary
        Add(new ILDParallelPlanarMeasLayer(silicon, air, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SITSenRear" )) ;  
      }
      
      
    }    
    
  }
  
  SetOwner();                   
}



void ILDSITKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  
  const gear::ZPlanarParameters& pSITDetMain = gearMgr.getSITParameters();
  const gear::ZPlanarLayerLayout& pSITLayerLayout = pSITDetMain.getZPlanarLayerLayout();
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  _nLayers = pSITLayerLayout.getNLayers(); 
  _SITgeo.resize(_nLayers);
  
  bool n_sensors_per_ladder_present = true;
  
  try {

    std::vector<int> v = pSITDetMain.getIntVals("n_sensors_per_ladder");

  } catch (gear::UnknownParameterException& e) {

    n_sensors_per_ladder_present = false;

  }

  double strip_angle_deg = 0.0;
  
  try {
    
    strip_angle_deg = pSITDetMain.getDoubleVal("strip_angle_deg");
    _isStripDetector = true;
    
  } catch (gear::UnknownParameterException& e) {
    
    _isStripDetector = false;
    
  }
  
  
  //SJA:FIXME: for now the support is taken as the same size the sensitive
  //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
  //           which if traversed in the reverse direction to the next boundary then the track would be propagated through carbon
  //           for a significant distance 
  //std::cout << "=============SIT strip angle: " << strip_angle_deg << "==============" << std::endl;
  for( int layer=0; layer < _nLayers; ++layer){
    
      
    _SITgeo[layer].nLadders = pSITLayerLayout.getNLadders(layer); 
    _SITgeo[layer].phi0 = pSITLayerLayout.getPhi0(layer); 
    _SITgeo[layer].dphi = 2*M_PI / _SITgeo[layer].nLadders; 
    _SITgeo[layer].senRMin = pSITLayerLayout.getSensitiveDistance(layer); 
    _SITgeo[layer].supRMin = pSITLayerLayout.getLadderDistance(layer); 
    _SITgeo[layer].length = pSITLayerLayout.getSensitiveLength(layer)*2.0; // note: gear for historical reasons uses the halflength 
    _SITgeo[layer].width = pSITLayerLayout.getSensitiveWidth(layer); 
    _SITgeo[layer].offset = pSITLayerLayout.getSensitiveOffset(layer); 
    _SITgeo[layer].senThickness = pSITLayerLayout.getSensitiveThickness(layer); 
    _SITgeo[layer].supThickness = pSITLayerLayout.getLadderThickness(layer); 

    if (n_sensors_per_ladder_present) {
      _SITgeo[layer].nSensorsPerLadder =   pSITDetMain.getIntVals("n_sensors_per_ladder")[layer];
    }
    else{
      _SITgeo[layer].nSensorsPerLadder = 1 ;
    }
    
    _SITgeo[layer].sensorLength = _SITgeo[layer].length / _SITgeo[layer].nSensorsPerLadder;
    
    
    if (_isStripDetector) {
      _SITgeo[layer].stripAngle = strip_angle_deg * M_PI/180 ;
    } else {
      _SITgeo[layer].stripAngle = 0.0 ;
    }
    //std::cout << _SITgeo[layer].nLadders << " " << _SITgeo[layer].phi0 << " "<< _SITgeo[layer].dphi << " " << _SITgeo[layer].senRMin << " " << _SITgeo[layer].supRMin << " "
    //          << _SITgeo[layer].length << " " << _SITgeo[layer].width << " " << _SITgeo[layer].offset << " " << _SITgeo[layer].senThickness << " " << _SITgeo[layer].supThickness << " "
    //          << _SITgeo[layer].nSensorsPerLadder << " " << _SITgeo[layer].sensorLength << std::endl;
    // streamlog_out(DEBUG0) << " layer  = " << layer << std::endl;
    // streamlog_out(DEBUG0) << " nSensorsPerLadder  = " << _SITgeo[layer].nSensorsPerLadder << std::endl;
    // streamlog_out(DEBUG0) << " sensorLength  = " << _SITgeo[layer].sensorLength << std::endl;
    // streamlog_out(DEBUG0) << " stripAngle  = " << _SITgeo[layer].stripAngle << std::endl;
    // streamlog_out(DEBUG0) << " _isStripDetector  = " << _isStripDetector << std::endl;
    
  }
  
}

void ILDSITKalDetector::setupGearGeom( IGeomSvc* geoSvc ){

  dd4hep::DetElement world = geoSvc->getDD4HepGeo();
  dd4hep::DetElement sit;
  const std::map<std::string, dd4hep::DetElement>& subs = world.children();
  for(std::map<std::string, dd4hep::DetElement>::const_iterator it=subs.begin();it!=subs.end();it++){
    if(it->first!="SIT") continue;
    sit = it->second;
  }
  dd4hep::rec::ZPlanarData* sitData = nullptr;
  try{
    sitData = sit.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    std::cout << e.what() << " " << sitData << std::endl;
    throw GaudiException(e.what(), "FATAL", StatusCode::FAILURE);
  }

  const dd4hep::Direction& field = geoSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
  _bZ = field.z()/dd4hep::tesla;

  std::vector<dd4hep::rec::ZPlanarData::LayerLayout>& sitlayers = sitData->layers;
  _nLayers = sitlayers.size();
  _SITgeo.resize(_nLayers);

  double strip_angle_deg = sitData->angleStrip/CLHEP::degree;
  if(strip_angle_deg!=0){
    _isStripDetector = true;
  }
  //std::cout << "=============SIT strip angle: " << strip_angle_deg << "==============" << std::endl;
  for( int layer=0; layer < _nLayers; ++layer){
    dd4hep::rec::ZPlanarData::LayerLayout& pSITLayerLayout = sitlayers[layer];
    
    _SITgeo[layer].nLadders = pSITLayerLayout.ladderNumber;
    _SITgeo[layer].phi0 = pSITLayerLayout.phi0;
    _SITgeo[layer].dphi = 2*M_PI / _SITgeo[layer].nLadders;
    _SITgeo[layer].senRMin = pSITLayerLayout.distanceSensitive*CLHEP::cm;
    _SITgeo[layer].supRMin = pSITLayerLayout.distanceSupport*CLHEP::cm;
    _SITgeo[layer].length = pSITLayerLayout.zHalfSensitive*2.0*CLHEP::cm; // note: gear for historical reasons uses the halflength                                                         
    _SITgeo[layer].width = pSITLayerLayout.widthSensitive*CLHEP::cm;
    _SITgeo[layer].offset = pSITLayerLayout.offsetSensitive*CLHEP::cm;
    _SITgeo[layer].senThickness = pSITLayerLayout.thicknessSensitive*CLHEP::cm;
    _SITgeo[layer].supThickness = pSITLayerLayout.thicknessSupport*CLHEP::cm;
    _SITgeo[layer].nSensorsPerLadder = pSITLayerLayout.sensorsPerLadder;
    _SITgeo[layer].sensorLength = _SITgeo[layer].length / _SITgeo[layer].nSensorsPerLadder;

    if (_isStripDetector) {
      _SITgeo[layer].stripAngle = strip_angle_deg * M_PI/180 ;
    } else {
      _SITgeo[layer].stripAngle = 0.0 ;
    }
    //std::cout << _SITgeo[layer].nLadders << " " << _SITgeo[layer].phi0 << " "<< _SITgeo[layer].dphi << " " << _SITgeo[layer].senRMin << " " << _SITgeo[layer].supRMin << " "
    //	      << _SITgeo[layer].length << " " << _SITgeo[layer].width << " " << _SITgeo[layer].offset << " " << _SITgeo[layer].senThickness << " " << _SITgeo[layer].supThickness << " "
    //	      << _SITgeo[layer].nSensorsPerLadder << " " << _SITgeo[layer].sensorLength << " " << pSITLayerLayout.lengthSensor << std::endl; 
  }
}
